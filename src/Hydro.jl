using HydroEOS, GridTypes, PyPlot
include("ParticleMesh.jl")

HydroAnalysis(gp, gd, p) = analyzePM(gp, gd, p)

function EvolveHydroEulerCosmology(gp, gd, p)
    @debug("Entered EvolveHydroEulerCosmology")
    if !isdefined(NoName,:cosmo)
        @warn("EvolveHydroEulerCosmology: No cosmology defined!")
    end
    if isdefined(NoName,:darkmatter) && darkmatter
        EvolveHydroEulerWithDMCosmology(gp, gd, p)
    else
        EvolveHydroEulerOnlyCosmology(gp, gd, p)
    end
end

function EvolveHydroEulerWithDMCosmology(gp, gd, p)
    @debug("Entered EvolveHydroEulerWithDMCosmology")
    a = 1.0
    ȧ = 0.0
    zinit = conf["InitialRedshift"]
    # Dark matter quantities
    g = gp[1]
    d = gd[1]
    rho = d.d["ρD"]         # DM density
    phi = d.d["Φ"]          # potential
    x = p["x"]
    v = p["v"]
    m = p["m"]
    acc   = zeros((3,g.dim...)) 
    pa = zeros(x)          # particle accelerations
    rt = rfft(phi) # initialize buffer for real-to complex transform

    # Baryon quantities
    Ng1 = g.dim[1]
    Ng2 = g.dim[2]
    Ng3 = g.dim[3]
    Neq = 5
    U = zeros(Neq, Ng1, Ng2, Ng3)
    U[1,:,:,:] = d.d["ρ"]
    U[2,:,:,:] = d.d["ρ"].*d.d["E"]
    U[3,:,:,:] = d.d["ρ"].*d.d["V₁"]
    U[4,:,:,:] = d.d["ρ"].*d.d["V₂"]
    U[5,:,:,:] = d.d["ρ"].*d.d["V₃"]

    dUF = zeros(Neq, Ng1, Ng2, Ng3)
    dUG = zeros(Neq, Ng1, Ng2, Ng3)
    dUH = zeros(Neq, Ng1, Ng2, Ng3)
    FF  = zeros(Neq, Ng1+1,Ng2,  Ng3)
    GG  = zeros(Neq, Ng1,  Ng2+1,Ng3)
    HH  = zeros(Neq, Ng1,  Ng2,  Ng3+1)

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    dt = InitialDt
    c = conf # conevenience
    
    c["CurrentTime"] = c["StartTime"]
    c["CurrentCycle"] = c["StartCycle"]    
    c["CurrentRedshift"] = z_from_t_internal(cosmo, c["CurrentTime"], zinit)
    ctr = 0
    tstep = (c["StopTime"]-c["CurrentTime"])/10

    while c["CurrentTime"] < c["StopTime"] && c["CurrentCycle"] < c["StopCycle"]
        if conf["CurrentTime"] >= tstep*ctr
            HydroAnalysis(gp, gd, p)
            ctr += 1
        end
        # gravitational evolution of DM + baryon
        make_periodic(x)
        deposit(rho,x,m,interpolation=ParticleDepositInterpolation)
        temprho = rho+d.d["ρ"] # add density from baryons
        rho_mean = mean(temprho)
        temprho -= rho_mean # the total sum of the density needs to be zero for a periodic signal
        compute_potential(phi, temprho, rt)
        compute_acceleration(acc, 100*phi)
        interpVecToPoints(pa, acc, x, interpolation=ParticleBackInterpolation)
        #println(pa[:,15])
        # LeapFrog step. Some codes combine the two drift steps 
        # note DM only here
        a = a_from_t_internal(cosmo, c["CurrentTime"] + dt/2/2, zinit, zwherea1=zinit)
        drift(x,v,dt/2/a)
        a = a_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zwherea1=zinit)
        ȧ = dadt_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zhwerea1=zinit)
        kick(v,pa,dt/a, ȧ)
        a = a_from_t_internal(cosmo, c["CurrentTime"] + 3.0*dt/2/2, zinit, zwherea1=zinit)
        drift(x,v,dt/2/a)

        # For baryon
        HLL_fluxes(FF, GG, HH, U, gd)
        update_dUF(dUF,FF) # don't do the 1./dx here.. : equ(3)
        update_dUG(dUG,GG)
        update_dUH(dUH,HH)
        a = a_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zwherea1=zinit)
        ȧ = dadt_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zhwerea1=zinit)
        update_UCosmology(U, dUF,dUG,dUH, dx, acc, dt/a, ȧ)

        primitive_from_U(gd, U)

        # update time
        c["CurrentTime"] += dt
        c["CurrentRedshift"] = z_from_t_internal(cosmo, c["CurrentTime"], zinit)
        dtFrac = DtFractionOfTime*c["CurrentTime"]
        dt = min(ComputeHydroTimeStep(gp, gd, p), CourantFactor*minimum(dx)/maximum(abs(p["v"])+1e-30), dtFrac, MaximumDt)
        if ((c["CurrentTime"]+dt) > c["StopTime"] )
            dt = (1. + 1e-6)*(c["StopTime"] - c["CurrentTime"])
            println("final dt = ", dt)
        end
        c["CurrentCycle"] += 1
    end
    @info("Evolved to time:", c["CurrentTime"])
    @info("Number of cycles:", c["CurrentCycle"])
    HydroAnalysis(gp, gd, p)

    nothing
end


function EvolveHydroEulerOnlyCosmology(gp, gd, p)
    @debug("Entered EvolveHydroEulerOnlyCosmology")
    # Baryon quantities
    g = gp[1]
    d = gd[1]
    Ng1 = g.dim[1]
    Ng2 = g.dim[2]
    Ng3 = g.dim[3]
    Neq = 5
    U = zeros(Neq, Ng1, Ng2, Ng3)
    U[1,:,:,:] = d.d["ρ"]
    U[2,:,:,:] = d.d["ρ"].*d.d["E"]
    U[3,:,:,:] = d.d["ρ"].*d.d["V₁"]
    U[4,:,:,:] = d.d["ρ"].*d.d["V₂"]
    U[5,:,:,:] = d.d["ρ"].*d.d["V₃"]

    dUF = zeros(Neq, Ng1, Ng2, Ng3)
    dUG = zeros(Neq, Ng1, Ng2, Ng3)
    dUH = zeros(Neq, Ng1, Ng2, Ng3)
    FF  = zeros(Neq, Ng1+1,Ng2,  Ng3)
    GG  = zeros(Neq, Ng1,  Ng2+1,Ng3)
    HH  = zeros(Neq, Ng1,  Ng2,  Ng3+1)

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    dt = InitialDt
    c = conf # conevenience
    zinit = c["InitialRedshift"]
    c["CurrentTime"] = c["StartTime"]
    c["CurrentCycle"] = c["StartCycle"]    
    c["CurrentRedshift"] = z_from_t_internal(cosmo, c["CurrentTime"], zinit)
    ctr = 0
    tstep = (c["StopTime"]-c["CurrentTime"])/10

    while c["CurrentTime"] < c["StopTime"] && c["CurrentCycle"] < c["StopCycle"]
        if conf["CurrentTime"] >= tstep*ctr
            @debug("Current cycle : ", string(conf["CurrentCycle"]))
            @debug("Current time : ", string(round(conf["CurrentTime"],5)))
            #HydroAnalysis(gp, gd, p)
            ctr += 1
        end
        
        HLL_fluxes(FF, GG, HH, U, gd)
     
        update_dUF(dUF,FF) # don't do the 1./dx here.. : equ(3)
        update_dUG(dUG,GG)
        update_dUH(dUH,HH)
        
        a = a_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zwherea1=zinit)
        ȧ = dadt_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zhwerea1=zinit)
        update_UCosmology(U, dUF,dUG,dUH, dx, acc, dt/a, ȧ)

        primitive_from_U(gd, U)

        # update time
        c["CurrentTime"] += dt
        c["CurrentRedshift"] = z_from_t_internal(cosmo, c["CurrentTime"], zinit)
        dtFrac = DtFractionOfTime*c["CurrentTime"]
        dt = min(ComputeHydroTimeStep(gp, gd, p), dtFrac, MaximumDt)
        if ((c["CurrentTime"]+dt) > c["StopTime"] )
            dt = (1. + 1e-6)*(c["StopTime"] - c["CurrentTime"])
            println("final dt = ", dt)
        end
        c["CurrentCycle"] += 1
    end
    @info("Evolved to time:", c["CurrentTime"])
    @info("Number of cycles:", c["CurrentCycle"])
    #HydroAnalysis(gp, gd, p)

    nothing
end

function EvolveHydroEuler(gp, gd, p)
    @debug("Entered EvolveHydroEuler")

    if isdefined(NoName,:cosmo)
        @warn("EvolveHydroEuler: You defined a cosmology. This Evolve routine will not use it.")
        @warn("EvolveHydroEuler: You may want to use EvolveHydroEulerCosmology")
    end
    if isdefined(NoName,:darkmatter) && darkmatter
        EvolveHydroEulerWithDM(gp, gd, p)
    else
        EvolveHydroEulerOnly(gp, gd, p)
    end
end

function EvolveHydroEulerWithDM(gp, gd, p)
    @debug("Entered EvolveHydroEulerWithDM")
    # Dark matter quantities
    g = gp[1]
    d = gd[1]
    rho = d.d["ρD"]         # DM density
    phi = d.d["Φ"]          # potential
    x = p["x"]
    v = p["v"]
    m = p["m"]
    acc   = zeros((3,g.dim...)) 
    pa = zeros(x)          # particle accelerations
    rt = rfft(phi) # initialize buffer for real-to complex transform

    # Baryon quantities
    Ng1 = g.dim[1]
    Ng2 = g.dim[2]
    Ng3 = g.dim[3]
    Neq = 5
    U = zeros(Neq, Ng1, Ng2, Ng3)
    U[1,:,:,:] = d.d["ρ"]
    U[2,:,:,:] = d.d["ρ"].*d.d["E"]
    U[3,:,:,:] = d.d["ρ"].*d.d["V₁"]
    U[4,:,:,:] = d.d["ρ"].*d.d["V₂"]
    U[5,:,:,:] = d.d["ρ"].*d.d["V₃"]

    dUF = zeros(Neq, Ng1, Ng2, Ng3)
    dUG = zeros(Neq, Ng1, Ng2, Ng3)
    dUH = zeros(Neq, Ng1, Ng2, Ng3)
    FF  = zeros(Neq, Ng1+1,Ng2,  Ng3)
    GG  = zeros(Neq, Ng1,  Ng2+1,Ng3)
    HH  = zeros(Neq, Ng1,  Ng2,  Ng3+1)

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    dt = InitialDt
    c = conf # conevenience
    
    c["CurrentTime"] = c["StartTime"]
    c["CurrentCycle"] = c["StartCycle"]    
    ctr = 0
    
    while c["CurrentTime"] < c["StopTime"] && c["CurrentCycle"] < c["StopCycle"]
        if conf["CurrentTime"] >= 0.1*ctr
            HydroAnalysis(gp, gd, p)
            ctr += 1
        end
        # gravitational evolution of DM + baryon
        make_periodic(x)
        deposit(rho,x,m,interpolation=ParticleDepositInterpolation)
        rho += d.d["ρ"] # add density from baryons
        rho_mean = mean(rho)
        rho -= rho_mean # the total sum of the density needs to be zero for a periodic signal
        compute_potential(phi, rho, rt)
        compute_acceleration(acc, phi)
        interpVecToPoints(pa, acc, x, interpolation=ParticleBackInterpolation)

        # LeapFrog step. Some codes combine the two drift steps 
        # note DM only here
        drift(x,v,dt/2)
        kick(v,pa,dt)
        drift(x,v,dt/2)

        # For baryon
        HLL_fluxes(FF, GG, HH, U, gd)
        update_dUF(dUF,FF) # don't do the 1./dx here.. : equ(3)
        update_dUG(dUG,GG)
        update_dUH(dUH,HH)
        update_U(U, dUF,dUG,dUH, dx, acc, dt)

        primitive_from_U(gd, U)

        # update time
        c["CurrentTime"] += dt
        dt = minimum([ComputeHydroTimeStep(gp, gd, p), 
                      1.0001*StopTime-c["CurrentTime"], MaximumDt])
        dt = minimum([dt, CourantFactor*minimum(dx)/maximum(abs(p["v"])+1e-30)])  # timestep from DM
        c["CurrentCycle"] += 1
    end
    @info("Evolved to time:", c["CurrentTime"])
    @info("Number of cycles:", c["CurrentCycle"])
    HydroAnalysis(gp, gd, p)

    nothing
end


function EvolveHydroEulerOnly(gp, gd, p)
    @debug("Entered EvolveHydroEulerOnly")
    # Baryon quantities
    g = gp[1]
    d = gd[1]
    Ng1 = g.dim[1]
    Ng2 = g.dim[2]
    Ng3 = g.dim[3]
    Neq = 5
    U = zeros(Neq, Ng1, Ng2, Ng3)
    U[1,:,:,:] = d.d["ρ"]
    U[2,:,:,:] = d.d["ρ"].*d.d["E"]
    U[3,:,:,:] = d.d["ρ"].*d.d["V₁"]
    U[4,:,:,:] = d.d["ρ"].*d.d["V₂"]
    U[5,:,:,:] = d.d["ρ"].*d.d["V₃"]

    dUF = zeros(Neq, Ng1, Ng2, Ng3)
    dUG = zeros(Neq, Ng1, Ng2, Ng3)
    dUH = zeros(Neq, Ng1, Ng2, Ng3)
    FF  = zeros(Neq, Ng1+1,Ng2,  Ng3)
    GG  = zeros(Neq, Ng1,  Ng2+1,Ng3)
    HH  = zeros(Neq, Ng1,  Ng2,  Ng3+1)

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    dt = InitialDt
    c = conf # conevenience
    
    c["CurrentTime"] = c["StartTime"]
    c["CurrentCycle"] = c["StartCycle"]    
    ctr = 0

    while c["CurrentTime"] < c["StopTime"] && c["CurrentCycle"] < c["StopCycle"]
        if conf["CurrentTime"] >= 0.1*ctr
            @debug("Current cycle : ", string(conf["CurrentCycle"]))
            @debug("Current time : ", string(round(conf["CurrentTime"],5)))
            #HydroAnalysis(gp, gd, p)
            ctr += 1
        end
        
        HLL_fluxes(FF, GG, HH, U, gd)
     
        update_dUF(dUF,FF) # don't do the 1./dx here.. : equ(3)
        update_dUG(dUG,GG)
        update_dUH(dUH,HH)
        
        update_U(U, dUF,dUG,dUH, dx, dt)

        primitive_from_U(gd, U)

        # update time
        c["CurrentTime"] += dt
        dt = minimum([ComputeHydroTimeStep(gp, gd, p), 
                      1.0001*StopTime-c["CurrentTime"], MaximumDt])
        c["CurrentCycle"] += 1
    end
    @info("Evolved to time:", c["CurrentTime"])
    @info("Number of cycles:", c["CurrentCycle"])
    #HydroAnalysis(gp, gd, p)

    nothing
end

function primitive_from_U(gd, U)
    gd[1].d["ρ"] = slice(U,1,:,:,:)
    gd[1].d["E"] = slice(U,2,:,:,:)./gd[1].d["ρ"]
    gd[1].d["V₁"] = slice(U,3,:,:,:)./gd[1].d["ρ"]
    gd[1].d["V₂"] = slice(U,4,:,:,:)./gd[1].d["ρ"]
    gd[1].d["V₃"] = slice(U,5,:,:,:)./gd[1].d["ρ"]
    nothing
end

function update_U(U, dUF,dUG,dUH, dx, dt)
    for i in eachindex(U)
        @inbounds U[i] += (dt/dx[1]) * dUF[i] + (dt/dx[2]) * dUG[i] + (dt/dx[3]) * dUH[i]
    end
    nothing
end

# update_U when there is a potential
function update_U(U, dUF,dUG,dUH, dx, acc, dt)
    for i in eachindex(U[1,:,:,:])
        @inbounds begin
        # store the current rho, v1, v2, v3 for using in updating
            rho = U[1,i]
            rhov1 = U[3,i]
            rhov2 = U[4,i]
            rhov3 = U[5,i]
            U[1,i] += (dt/dx[1])*dUF[1,i] + (dt/dx[2])*dUG[1,i] + (dt/dx[3])*dUH[1,i]
            U[2,i] += (dt/dx[1])*dUF[2,i] + (dt/dx[2])*dUG[2,i] + (dt/dx[3])*dUH[2,i] 
                        + dt/2*((U[3,i]+rhov1)*acc[1,i]+(U[4,i]+rhov2)*acc[2,i]+(U[5,i]+rhov3)*acc[3,i])
            U[3,i] += (dt/dx[1])*dUF[3,i] + (dt/dx[2])*dUG[3,i] + (dt/dx[3])*dUH[3,i] + dt/2*(U[1,i]+rho)*acc[1,i]
            U[4,i] += (dt/dx[1])*dUF[4,i] + (dt/dx[2])*dUG[4,i] + (dt/dx[3])*dUH[4,i] + dt/2*(U[1,i]+rho)*acc[2,i]
            U[5,i] += (dt/dx[1])*dUF[5,i] + (dt/dx[2])*dUG[5,i] + (dt/dx[3])*dUH[5,i] + dt/2*(U[1,i]+rho)*acc[3,i]
        end
    end
    nothing
end

# update_U when there is a potential + cosmology
function update_UCosmology(U, dUF,dUG,dUH, dx, acc, dt, dadt)
    #semi implicit integration
    coef = 0.5*dadt*dt
    coef1 = 1.0 - coef
    coef2 = 1.0 / (1.0 + coef)
    coefE = dadt*dt
    coef1E = 1.0 - coefE
    coef2E = 1.0 / (1.0 + coefE)
    for i in eachindex(U[1,:,:,:])
        @inbounds begin
        # store the current rho, v1, v2, v3 for using in updating
            rho = U[1,i]
            rhov1 = U[3,i]
            rhov2 = U[4,i]
            rhov3 = U[5,i]
            U[1,i] += (dt/dx[1])*dUF[1,i] + (dt/dx[2])*dUG[1,i] + (dt/dx[3])*dUH[1,i]
            U[2,i] = (coef1E*U[2,i] + (dt/dx[1])*dUF[2,i] + (dt/dx[2])*dUG[2,i] + (dt/dx[3])*dUH[2,i]
                        + dt/2*((U[3,i]+rhov1)*acc[1,i]+(U[4,i]+rhov2)*acc[2,i]+(U[5,i]+rhov3)*acc[3,i]))*coef2E
            U[3,i] = (coef1*U[3,i] + (dt/dx[1])*dUF[3,i] + (dt/dx[2])*dUG[3,i] + (dt/dx[3])*dUH[3,i] + dt/2*(U[1,i]+rho)*acc[1,i])*coef2
            U[4,i] = (coef1*U[4,i] + (dt/dx[1])*dUF[4,i] + (dt/dx[2])*dUG[4,i] + (dt/dx[3])*dUH[4,i] + dt/2*(U[1,i]+rho)*acc[2,i])*coef2
            U[5,i] = (coef1*U[5,i] + (dt/dx[1])*dUF[5,i] + (dt/dx[2])*dUG[5,i] + (dt/dx[3])*dUH[5,i] + dt/2*(U[1,i]+rho)*acc[3,i])*coef2
        end
    end
    nothing
end

# update_U with no potential, but cosmology
function update_UCosmology(U, dUF,dUG,dUH, dx, dt, dadt)
    #semi implicit integration
    coef = 0.5*dadt*dt
    coef1 = 1.0 - coef
    coef2 = 1.0 / (1.0 + coef)
    coefE = dadt*dt
    coef1E = 1.0 - coefE
    coef2E = 1.0 / (1.0 + coefE)
    for i in eachindex(U[1,:,:,:])
        @inbounds begin
        # store the current rho, v1, v2, v3 for using in updating
            rho = U[1,i]
            rhov1 = U[3,i]
            rhov2 = U[4,i]
            rhov3 = U[5,i]
            U[1,i] += (dt/dx[1])*dUF[1,i] + (dt/dx[2])*dUG[1,i] + (dt/dx[3])*dUH[1,i]
            U[2,i] = (coef1E*U[2,i] + (dt/dx[1])*dUF[2,i] + (dt/dx[2])*dUG[2,i] + (dt/dx[3])*dUH[2,i])*coef2E
            U[3,i] = (coef1*U[3,i] + (dt/dx[1])*dUF[3,i] + (dt/dx[2])*dUG[3,i] + (dt/dx[3])*dUH[3,i])*coef2
            U[4,i] = (coef1*U[4,i] + (dt/dx[1])*dUF[4,i] + (dt/dx[2])*dUG[4,i] + (dt/dx[3])*dUH[4,i])*coef2
            U[5,i] = (coef1*U[5,i] + (dt/dx[1])*dUF[5,i] + (dt/dx[2])*dUG[5,i] + (dt/dx[3])*dUH[5,i])*coef2
        end
    end
    nothing
end


function update_dUF(dUF,FF)
    for i=1:size(dUF,2), n=1:size(dUF,1)
        @inbounds dUF[n,i,:,:] = -(FF[n,i+1,:,:] - FF[n,i,:,:]) # so this is -(F_{j+1/2}-F_{j-1/2})
    end
    nothing
end

function update_dUG(dUG,GG)
    for j=1:size(dUG,3), n=1:size(dUG,1)
        @inbounds dUG[n,:,j,:] = -(GG[n,:,j+1,:] - GG[n,:,j,:]) 
    end
    nothing
end

function update_dUH(dUH,HH)
    for k=1:size(dUH,4), n=1:size(dUH,1)
        @inbounds dUH[n,:,:,k] = -(HH[n,:,:,k+1] - HH[n,:,:,k])
    end
    nothing
end


function HLL_fluxes(FF, GG, HH, U, gd)
    d = gd[1]
    Ngl = size(U)[2:end]
    Neq = size(U)[1]
    UL_F = zeros(Neq, Ngl[1]+1, Ngl[2],Ngl[3])
    UR_F = zeros(UL_F)
    FL = zeros(UL_F)
    FR = zeros(UL_F)
    UL_G = zeros(Neq, Ngl[1], Ngl[2]+1,Ngl[3])
    UR_G = zeros(UL_G)
    GL = zeros(UL_G)
    GR = zeros(UL_G)
    UL_H = zeros(Neq, Ngl[1], Ngl[2],Ngl[3]+1)
    UR_H = zeros(UL_H)
    HL = zeros(UL_H)
    HR = zeros(UL_H)


    # create left and right states
    LR_states(UL_F,UR_F,UL_G,UR_G,UL_H,UR_H, U)


    # calculate P and cs
    v2 = zeros(1,Ngl...) ## note the weird indexing here.... for easier algebra with U 
    v2[1,:,:,:] = d.d["V₁"].^2 + d.d["V₂"].^2 + d.d["V₃"].^2
    if DualEnergyFormalism
        ε = d.d["ε"]
    else
        ε = zeros(1,Ngl...)
        ε[1,:,:,:] = U[2,:,:,:]./U[1,:,:,:] - 0.5*v2[1,:,:,:]
    end
    P = zeros(1,Ngl...)
    cs = zeros(1,Ngl...)
    for i=1:Ngl[1], j=1:Ngl[2], k=1:Ngl[3]
        P[1,i,j,k] = P_from_ε(ε[1,i,j,k], U[1,i,j,k], eos)
        cs[1,i,j,k] = cs_from_ε(ε[1,i,j,k], eos)
    end

    # create left and right fluxes
    LR_fluxes(FL,FR, GL,GR, HL,HR, U, gd, P)
    

    # create HLL fluxes, FF, GG, HH
    for n=1:size(FF,1), j=1:size(FF,2)
        @inbounds begin
            jm1 = mod1(j-1,Ngl[1])
            jj = mod1(j,Ngl[1])
            lam_l_p = sqrt(v2[1,jm1,:,:]) + cs[1,jm1,:,:]
            lam_l_m = sqrt(v2[1,jm1,:,:]) - cs[1,jm1,:,:]
            lam_r_p = sqrt(v2[1,jj,:,:]) + cs[1,jj,:,:]
            lam_r_m = sqrt(v2[1,jj,:,:]) - cs[1,jj,:,:]
            ap = zeros(lam_l_p)
            am = zeros(lam_l_p)
            for k in eachindex(lam_l_p)
                ap[k] = maximum([0.,lam_l_p[k],lam_r_p[k]])
                am[k] = maximum([0.,-lam_l_m[k],-lam_r_m[k]])
            end
            FF[n,j,:,:] = ((ap.*FL[n,j,:,:]+am.*FR[n,j,:,:]-ap.*am.*(UR_F[n,j,:,:]-UL_F[n,j,:,:])))./(ap+am)
        end
    end
    for n=1:size(GG,1), j=1:size(GG,3)
        @inbounds begin
            jm1 = mod1(j-1,Ngl[2])
            jj = mod1(j,Ngl[2])
            lam_l_p = sqrt(v2[1,:,jm1,:]) + cs[1,:,jm1,:]
            lam_l_m = sqrt(v2[1,:,jm1,:]) - cs[1,:,jm1,:]
            lam_r_p = sqrt(v2[1,:,jj,:]) + cs[1,:,jj,:]
            lam_r_m = sqrt(v2[1,:,jj,:]) - cs[1,:,jj,:]
            ap = zeros(lam_l_p)
            am = zeros(lam_l_p)
            for k in eachindex(lam_l_p)
                ap[k] = maximum([0.,lam_l_p[k],lam_r_p[k]])
                am[k] = maximum([0.,-lam_l_m[k],-lam_r_m[k]])
            end
            GG[n,:,j,:] = ((ap.*GL[n,:,j,:]+am.*GR[n,:,j,:]-ap.*am.*(UR_G[n,:,j,:]-UL_G[n,:,j,:])))./(ap+am)
        end
    end
    for n=1:size(HH,1), j=1:size(HH,4)
        @inbounds begin
            jm1 = mod1(j-1,Ngl[3])
            jj = mod1(j,Ngl[3])
            lam_l_p = sqrt(v2[1,:,:,jm1]) + cs[1,:,:,jm1]
            lam_l_m = sqrt(v2[1,:,:,jm1]) - cs[1,:,:,jm1]
            lam_r_p = sqrt(v2[1,:,:,jj])  + cs[1,:,:,jj]
            lam_r_m = sqrt(v2[1,:,:,jj])  - cs[1,:,:,jj]
            ap = zeros(lam_l_p)
            am = zeros(lam_l_p)
            for k in eachindex(lam_l_p)
                ap[k] = maximum([0.,lam_l_p[k],lam_r_p[k]])
                am[k] = maximum([0.,-lam_l_m[k],-lam_r_m[k]])
            end
            HH[n,:,:,j] = ((ap.*HL[n,:,:,j]+am.*HR[n,:,:,j]-ap.*am.*(UR_H[n,:,:,j]-UL_H[n,:,:,j])))./(ap+am)
        end
    end

    nothing
end

function LR_states(UL_F,UR_F,UL_G,UR_G,UL_H,UR_H, U)
    # piecewise constant 
    Ngl = size(U)[2:end]
    Neq = size(U)[1]
    for n=1:Neq
        for i=1:Ngl[1]+1, j=1:Ngl[2], k=1:Ngl[3]
            @inbounds UL_F[n,i,j,k] = U[n,mod1(i-1,end), j, k] # left is -1 : UL stores information of [end,1,2,...,end]
            @inbounds UR_F[n,i,j,k] = U[n,mod1(i,end), j, k]  # right is here : UR stores information of [1,2,...,end,1]
        end
    end
    for n=1:size(U,1)
        for i=1:Ngl[1], j=1:Ngl[2]+1, k=1:Ngl[3]
            @inbounds UL_G[n,i,j,k] = U[n,i,mod1(j-1,end), k] # left is -1 : UL stores information of [end,1,2,...,end]
            @inbounds UR_G[n,i,j,k] = U[n,i,mod1(j,end), k]  # right is here : UR stores information of [1,2,...,end,1]
        end
    end
    for n=1:size(U,1)
        for i=1:Ngl[1], j=1:Ngl[2], k=1:Ngl[3]+1
            @inbounds UL_H[n,i,j,k] = U[n,i,j,mod1(k-1,end)] # left is -1 : UL stores information of [end,1,2,...,end]
            @inbounds UR_H[n,i,j,k] = U[n,i,j,mod1(k,end)]  # right is here : UR stores information of [1,2,...,end,1]
        end
    end
    nothing
end

function LR_fluxes(FL,FR, GL,GR, HL,HR, U, gd, P)
    # piecewise constant 
    d = gd[1]
    Ngl = size(U)[2:end]
    Neq = size(U)[1]
    V₁ = zeros(1,Ngl...)
    V₂ = zeros(1,Ngl...)
    V₃ = zeros(1,Ngl...)
    V₁[1,:,:,:] = d.d["V₁"]
    V₂[1,:,:,:] = d.d["V₂"]
    V₃[1,:,:,:] = d.d["V₃"]

    
    for i=1:size(FL,2)
        @inbounds begin
            FL[1,i,:,:] = U[3,mod1(i-1,end),:,:]
            FR[1,i,:,:] = U[3,mod1(i,end),:,:]
            FL[2,i,:,:] = (U[2,mod1(i-1,end),:,:]+P[1,mod1(i-1,end),:,:]).*V₁[1,mod1(i-1,end),:,:]
            FR[2,i,:,:] = (U[2,mod1(i,end),:,:]+P[1,mod1(i,end),:,:]).*V₁[1,mod1(i,end),:,:]
            FL[3,i,:,:] = U[3,mod1(i-1,end),:,:].*V₁[1,mod1(i-1,end),:,:]+P[1,mod1(i-1,end),:,:]
            FR[3,i,:,:] = U[3,mod1(i,end),:,:].*V₁[1,mod1(i,end),:,:]+P[1,mod1(i,end),:,:]
            FL[4,i,:,:] = U[3,mod1(i-1,end),:,:].*V₂[1,mod1(i-1,end),:,:]
            FR[4,i,:,:] = U[3,mod1(i,end),:,:].*V₂[1,mod1(i,end),:,:]
            FL[5,i,:,:] = U[3,mod1(i-1,end),:,:].*V₃[1,mod1(i-1,end),:,:]
            FR[5,i,:,:] = U[3,mod1(i,end),:,:].*V₃[1,mod1(i,end),:,:]
        end
    end
    for i=1:size(GL,3)
        @inbounds begin
            GL[1,:,i,:] = U[4,:,mod1(i-1,end),:]
            GR[1,:,i,:] = U[4,:,mod1(i,end),:]
            GL[2,:,i,:] = (U[2,:,mod1(i-1,end),:]+P[1,:,mod1(i-1,end),:]).*V₂[1,:,mod1(i-1,end),:]
            GR[2,:,i,:] = (U[2,:,mod1(i,end),:]+P[1,:,mod1(i,end),:]).*V₂[1,:,mod1(i,end),:]
            GL[3,:,i,:] = U[4,:,mod1(i-1,end),:].*V₁[1,:,mod1(i-1,end),:]
            GR[3,:,i,:] = U[4,:,mod1(i,end),:].*V₁[1,:,mod1(i,end),:]
            GL[4,:,i,:] = U[4,:,mod1(i-1,end),:].*V₂[1,:,mod1(i-1,end),:]+P[1,:,mod1(i-1,end),:]
            GR[4,:,i,:] = U[4,:,mod1(i,end),:].*V₂[1,:,mod1(i,end),:]+P[1,:,mod1(i,end),:]
            GL[5,:,i,:] = U[4,:,mod1(i-1,end),:].*V₃[1,:,mod1(i-1,end),:]
            GR[5,:,i,:] = U[4,:,mod1(i,end),:].*V₃[1,:,mod1(i,end),:]
        end
    end
    for i=1:size(HL,4)
        @inbounds begin
            HL[1,:,:,i] = U[5,:,:,mod1(i-1,end)]
            HR[1,:,:,i] = U[5,:,:,mod1(i,end)]
            HL[2,:,:,i] = (U[2,:,:,mod1(i-1,end)]+P[1,:,:,mod1(i-1,end)]).*V₃[1,:,:,mod1(i-1,end)]
            HR[2,:,:,i] = (U[2,:,:,mod1(i,end)]+P[1,:,:,mod1(i,end)]).*V₃[1,:,:,mod1(i,end)]
            HL[3,:,:,i] = U[5,:,:,mod1(i-1,end)].*V₁[1,:,:,mod1(i-1,end)]
            HR[3,:,:,i] = U[5,:,:,mod1(i,end)].*V₁[1,:,:,mod1(i,end)]
            HL[4,:,:,i] = U[5,:,:,mod1(i-1,end)].*V₂[1,:,:,mod1(i-1,end)]
            HR[4,:,:,i] = U[5,:,:,mod1(i,end)].*V₂[1,:,:,mod1(i,end)]
            HL[5,:,:,i] = U[5,:,:,mod1(i-1,end)].*V₃[1,:,:,mod1(i-1,end)]+P[1,:,:,mod1(i-1,end)]
            HR[5,:,:,i] = U[5,:,:,mod1(i,end)].*V₃[1,:,:,mod1(i,end)]+P[1,:,:,mod1(i,end)]
        end
    end
    nothing
end

function ComputeHydroTimeStep(gp, gd, p)
    # calculate P and cs 
    d = gd[1]
    g = gp[1]
    Ngl = size(d.d["ρ"])
    v2 = d.d["V₁"].^2 + d.d["V₂"].^2 + d.d["V₃"].^2
    if DualEnergyFormalism
        ε = d.d["ε"]
    else
        ε = d.d["E"] - 0.5*v2
    end
    cs = zeros(Ngl...)
    for k=1:Ngl[3], j=1:Ngl[2], i=1:Ngl[1]
        cs[i,j,k] = cs_from_ε(ε[i,j,k], eos)
    end

    max_cs = maximum(cs)
    max_v = maximum(sqrt(v2))
    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    return CourantFactor * minimum(dx)/(max_cs+max_v)
end

function InitializeHydroSimulation(gp, gd, p)
    @debug("InitializeHydroSimulation: start.")
    global const rank = sum([(i > 1) for i in TopgridDimensions]) # infer dimensionality 
    global const eos = EOS(γ)
    if verbose
        println("Dimensions according to TopgridDimensions: ", rank)
    end
    
    # Create TopGrid
    dim =  ntuple((i) -> TopgridDimensions[i], length(TopgridDimensions))
    g = gp[1] = GridPatch(DomainLeftEdge,DomainRightEdge,1,1, dim) 

    gd[1] = GridData(Dict())
    d = gd[1]
    
    # allocate our arrays for the hydro quantities
    d.d["ρ"]  = ones(g.dim)
    d.d["E"]  = ones(g.dim)
    if DualEnergyFormalism
        d.d["ε"]  = ones(g.dim)
    end
    d.d["V₁"] = zeros(g.dim)
    d.d["V₂"] = zeros(g.dim)
    d.d["V₃"] = zeros(g.dim)

    #dictionary_to_variable_names(d.d) # shortcuts for the variable names

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    if darkmatter
        @debug("Initialize DM part: start.")
        Npart = prod(ParticleDimensions)

        # allocate our arrays for the grid quantities
        d.d["ρD"]  = ones(gp[1].dim) # mass density 
        d.d["Φ"]   = ones(gp[1].dim) # gravitational potential

        # allocate particles
        p["x"] = zeros(rank, Npart)
        p["v"] = zeros(rank, Npart)
        p["m"] = 1. * length(d.d["ρD"])/Npart # unity for particle mass.
    end    

    # define cosmology if it is specified in parameter file
    if haskey(conf, "Cosmology")
        eval(parse(string("global const cosmo=",conf["Cosmology"])))
        if haskey(conf, "InitialRedshift") # set start times from redshift
            zstart = conf["InitialRedshift"]
            u = get_units(cosmo, zstart, zinit=zstart)
            ε0 = 1.5*1.38064852e-16*u.Temp            
            tstart = age_gyr(cosmo, zstart)*(gyr_in_s/u.T)
            zend = conf["FinalRedshift"]
            tend = age_gyr(cosmo, zend)*(gyr_in_s/u.T)
            conf["CurrentTime"] = conf["StartTime"] = tstart
            conf["StopTime"] = tend
            @info("Start Redshift:", zstart, " -> t_start = ",tstart)
            @info("Final Redshift:", zend, " -> t_final = ",tend)

            conf["OmegaCDM"] = cosmo.Ω_m
            if haskey(conf,"OmegaBaryon")
                conf["OmegaCDM"] = cosmo.Ω_m - conf["OmegaBaryon"]
                d.d["ρ"] .*= conf["OmegaBaryon"] 
                d.d["E"] .*= 1. #conf["OmegaBaryon"] # / u.ρ *ε0
            else
                d.d["ρ"] .*= 0.2
                d.d["E"] .*= 0.2 # / u.ρ *ε0
            end         
            p["m"] .*= conf["OmegaCDM"] # particles
        else
            @warn("Specified Cosmology but gave no InitialRedshift! ")
            @warn("Proceed with caution ...")
        end
    end

     # call Initializer given as ProblemDefinition in the .conf file
    callsig = string("initialvalues=",conf["ProblemDefinition"])
    eval(parse(callsig)) # call the correct routine
    
    initialvalues(gp, gd, p)
    nothing
end

function SodShocktube_init()

    # Use this gamma law hydro equation of state
    # global const eos = EOS(γ) -- moved into function InitializeHydroSimulation(gp, gd, p)

    # Some basic error checking
    if maximum(DomainRightEdge) > 1 || minimum(DomainLeftEdge) > 1
        error("Unexpected DomainSize")
    end
    if eos.γ != 1.4
        println("Warning: Sod test traditionaly is with γ=1.4.")
        println("         We have γ=",eos.γ)
    end
    1 # successfully initialized
end

function SodShocktube(gp, gd, p)
    if !isdefined(:SodInitialized)
        global const SodInitialized = SodShocktube_init()
    end

    g = gp[1]
    d = gd[1]
    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]

    ρ  = d.d["ρ"] 
    E  = d.d["E"] 
    V₁ = d.d["V₁"]
    V₂ = d.d["V₂"]
    V₃ = d.d["V₃"]

    Pl = 1
    Pr = 0.1795
    ρl = 0.25
    ρr = 1
    v1 = 0
    v2 = 0
    v3 = 0

    for k in 1:size(ρ,3)
        for j in 1:size(ρ,2)
            for i in 1:size(ρ,1)
                x = g.LE + 0.5*dx .* [i, j, k] # position
                r = [i,j,k]

                if x[1] < 0.5
                    ρ[r...] = ρl
                    E[r...] = ε_from_P(Pl, ρ[r...], eos) + 0.5*(v1^2+v2^2+v3^2) 
                else
                    ρ[r...] = ρr
                    E[r...] = ε_from_P(Pr, ρ[r...], eos) + 0.5*(v1^2+v2^2+v3^2) 
                end
                V₁[r...] = v1
                V₂[r...] = v2
                V₃[r...] = v3

                if DualEnergyFormalism
                    ε[r...] = E[r...] - (0.5 *
                                   (V₁[r...].^2+V₂[r...].^2+V₃[r...].^2))
                end
            end
        end
    end
end

function UniformBaryons(gp,gd,p)
    @debug("UniformBaryons: start.")
    g = gp[1]
    d = gd[1]
    
    P = 1. # temporary
    rho = 1.#/ prod(g.dim)

    ρ  = d.d["ρ"] = rho*ones(g.dim)
    E  = d.d["E"] = ones(g.dim)
    if DualEnergyFormalism
        d.d["ε"]  = ones(g.dim)
        ε = d.d["ε"]
    end
    V₁ = d.d["V₁"] = zeros(g.dim)
    V₂ = d.d["V₂"] = zeros(g.dim)
    V₃ = d.d["V₃"] = zeros(g.dim)

    

    for k in 1:size(ρ,3)
        for j in 1:size(ρ,2)
            for i in 1:size(ρ,1)
                r = [i,j,k]
                E[r...] = ε_from_P(P, ρ[r...], eos) + (0.5 *(V₁[r...].^2+V₂[r...].^2+V₃[r...].^2))
                if DualEnergyFormalism
                    ε[r...] = E[r...] - (0.5 *(V₁[r...].^2+V₂[r...].^2+V₃[r...].^2))
                end
            end
        end
    end
end
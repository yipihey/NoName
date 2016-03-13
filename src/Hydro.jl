using HydroEOS, GridTypes

function EvolveHydroEuler(gp, gd, p)
    @debug("Entered EvolveHydroEuler")

    dt = InitialDt
    CurrentTime = StartTime
    CurrentCycle = StartCycle
    while (CurrentTime < StopTime) && (CurrentCycle < StopCycle)
        SolveHydroEquations!(gp, gd, p, dt) # keep p for doing tracer particles
        CurrentTime += dt
        CurrentCycle += 1
        dt = minimum([ComputeHydroTimeStep(gp, gd, p),
                      1.0001*StopTime-CurrentTime])
#        HydroAnalysis(gp, gd, p)
    end
    @info("Evolved to time:", CurrentTime)

    nothing
end

function SolveHydroEquations!(gp, gd, p, dt)

    d = gd[1]
    g = gp[1]
    Ngl = size(d.d["ρ"])
    Neq = 5
    U = zeros(Neq, Ngl...)
    U[1,:,:,:] = d.d["ρ"]
    U[2,:,:,:] = d.d["ρ"].*d.d["E"]
    U[3,:,:,:] = d.d["ρ"].*d.d["V₁"]
    U[4,:,:,:] = d.d["ρ"].*d.d["V₂"]
    U[5,:,:,:] = d.d["ρ"].*d.d["V₃"]

    dUF = zeros(Neq, Ngl...)
    dUG = zeros(Neq, Ngl...)
    dUH = zeros(Neq, Ngl...)
    FF = zeros(Neq, Ngl[1]+1,Ngl[2],Ngl[3])
    GG = zeros(Neq, Ngl[1],Ngl[2]+1,Ngl[3])
    HH = zeros(Neq, Ngl[1],Ngl[2],Ngl[3]+1)

    HLL_fluxes(FF, GG, HH, U, gd)

    update_dUF(dUF,FF) # don't do the 1./dx here.. : equ(3)
    update_dUG(dUG,GG)
    update_dUH(dUH,HH)

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    update_U(U, dUF,dUG,dUH, dx, dt)

    primitive_from_U(gp, gd, p, U)
    nothing
end

function primitive_from_U(gp, gd, p, U)
    gd[1].d["ρ"] = slice(U,1,:,:,:)
    gd[1].d["E"] = slice(U,2,:,:,:)./gd[1].d["ρ"]
    gd[1].d["V₁"] = slice(U,3,:,:,:)./gd[1].d["ρ"]
    gd[1].d["V₂"] = slice(U,4,:,:,:)./gd[1].d["ρ"]
    gd[1].d["V₃"] = slice(U,5,:,:,:)./gd[1].d["ρ"]
    nothing
end

function update_U(U, dUF,dUG,dUH, dx, dt)
    for i in eachindex(U)
        U[i] += (dt/dx[1]) * dUF[i] + (dt/dx[2]) * dUG[i] + (dt/dx[3]) * dUH[i]
    end
    nothing
end


function update_dUF(dUF,FF)
    for i=1:size(dUF,2), n=1:size(dUF,1)
        dUF[n,i,:,:] = -(FF[n,i+1,:,:] - FF[n,i,:,:]) # so this is -(F_{j+1/2}-F_{j-1/2})
    end
    nothing
end

function update_dUG(dUG,GG)
    for j=1:size(dUG,3), n=1:size(dUG,1)
        dUG[n,:,j,:] = -(GG[n,:,j+1,:] - GG[n,:,j,:])
    end
    nothing
end

function update_dUH(dUH,HH)
    for k=1:size(dUH,4), n=1:size(dUH,1)
        dUH[n,:,:,k] = -(HH[n,:,:,k+1] - HH[n,:,:,k])
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
    for n=1:size(GG,1), j=1:size(GG,3)
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
    for n=1:size(HH,1), j=1:size(HH,4)
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

    nothing
end

function LR_states(UL_F,UR_F,UL_G,UR_G,UL_H,UR_H, U)
    # piecewise constant
    Ngl = size(U)[2:end]
    Neq = size(U)[1]
    for n=1:Neq
        for i=1:Ngl[1]+1, j=1:Ngl[2], k=1:Ngl[3]
            UL_F[n,i,j,k] = U[n,mod1(i-1,end), j, k] # left is -1 : UL stores information of [end,1,2,...,end]
            UR_F[n,i,j,k] = U[n,mod1(i,end), j, k]  # right is here : UR stores information of [1,2,...,end,1]
        end
    end
    for n=1:size(U,1)
        for i=1:Ngl[1], j=1:Ngl[2]+1, k=1:Ngl[3]
            UL_G[n,i,j,k] = U[n,i,mod1(j-1,end), k] # left is -1 : UL stores information of [end,1,2,...,end]
            UR_G[n,i,j,k] = U[n,i,mod1(j,end), k]  # right is here : UR stores information of [1,2,...,end,1]
        end
    end
    for n=1:size(U,1)
        for i=1:Ngl[1], j=1:Ngl[2], k=1:Ngl[3]+1
            UL_H[n,i,j,k] = U[n,i,j,mod1(k-1,end)] # left is -1 : UL stores information of [end,1,2,...,end]
            UR_H[n,i,j,k] = U[n,i,j,mod1(k,end)]  # right is here : UR stores information of [1,2,...,end,1]
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
    for i=1:size(GL,3)
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
    for i=1:size(HL,4)
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
    P = zeros(Ngl...)
    cs = zeros(Ngl...)
    for i=1:Ngl[1], j=1:Ngl[2], k=1:Ngl[3]
        P[i,j,k] = P_from_ε(ε[i,j,k], d.d["ρ"][i,j,k], eos)
        cs[i,j,k] = cs_from_ε(ε[i,j,k], eos)
    end

    max_cs = maximum(cs)
    max_v = maximum(sqrt(v2))
    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]
    return CourantFactor * minimum(dx)/(max_cs+max_v)
end

function InitializeHydroSimulation(gp, gd, p)

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

    ρ  = d.d["ρ"]
    E  = d.d["E"]
    V₁ = d.d["V₁"]
    V₂ = d.d["V₂"]
    V₃ = d.d["V₃"]

    #dictionary_to_variable_names(d.d) # shortcuts for the variable names

    dx = [(g.RE[i]-g.LE[i])/g.dim[i] for i in 1:3]

     # call Initializer given as ProblemDefinition in the .conf file
    callsig = string("initialvalues=",conf["ProblemDefinition"])
    eval(parse(callsig)) # call the correct routine

    for k in 1:size(ρ,3)
        for j in 1:size(ρ,2)
            for i in 1:size(ρ,1)
                x = g.LE + 0.5*dx .* [i, j, k] # position
                r = [i,j,k]
                (ρ[r], E[r], V₁[r], V₂[r], V₃[r]) =
                    initialvalues(x)
                if DualEnergyFormalism
                    ε[r] = E[r] - (0.5 *
                                   (V₁[r].^2+V₂[r].^2+V₃[r].^2))
                end
            end
        end
    end

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

function SodShocktube(x)
    if !isdefined(:SodInitialized)
        global const SodInitialized = SodShocktube_init()
    end

    Pl = 1
    Pr = 0.1795
    ρl = 0.25
    ρr = 1
    V₁ = 0
    V₂ = 0
    V₃ = 0
    if x[1] < 0.5
        ρ = ρl
        E = ε_from_P(Pl, ρ, eos) + 0.5*(V₁^2+V₂^2+V₃^2)
    else
        ρ = ρr
        E = ε_from_P(Pr, ρ, eos) + 0.5*(V₁^2+V₂^2+V₃^2)
    end
    ρ, E, V₁, V₂, V₃
end

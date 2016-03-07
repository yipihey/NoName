using PyPlot
include("InputPowerSpectra.jl")

function analyzePM(gp, gd, p)

    if !IOtools.readyForOutput(conf)
        return nothing
    end
    c = conf

    if !isdir(OutputDirectory)
        mkdir(OutputDirectory)
    end

    @info("Output:", c["CurrentOutputNumber"], ": CurrentTime:", c["CurrentTime"])
    if isdefined(:cosmo)
        @info("CurrentRedshift:", c["CurrentRedshift"])
    end
    
    on = c["CurrentOutputNumber"] # Since we are changing it cannot use shortcut
    s = @sprintf("%5.5i", on)
    newD = OutputDirectory
    if OutputSeparateDirectories
        newD = string(OutputDirectory,"/",OutputDirPrefix,s)
        if !isdir(newD)
            mkdir(newD)
        end
        @info(newD, " directory opened to write output")
    end

    # Write out parameters for easy reading
    fname = string(newD,"/",s,"_parameters.info")
    IOtools.outputGlobalParametersText(fname,conf)
    c["CurrentOutputNumber"] += 1

    # Write out grid and particle data
    fname = string(newD,"/",s,".h5")
    IOtools.grid_output(fname, gp, gd,overwrite=true)
    IOtools.particle_output(fname, p)
    
    
    # saving power spectrum
    fname = string(newD,"/",s,"_PS.png")
    karray, power = powerSpectrum(gp, gd, p)
    loglog(karray, power)
    xlabel("k")
    ylabel("P(k)")
    title("Power spectrum")
    savefig(fname)
    close() 
end

function InitializeParticleSimulation(gp, gd, p)
    @debug("InitializeDarkMatterSimulation: start.")
    global const rank = sum([(i > 1) for i in TopgridDimensions]) # infer dimensionality 
    @info("Dimensions according to TopgridDimensions: ", rank)
    
    Npart = prod(ParticleDimensions)

    # Create TopGrid
    dim =  (TopgridDimensions...)
    gp[1] = GridPatch(DomainLeftEdge,DomainRightEdge,1,1, dim) 

    gd[1] = GridData(Dict())
    d = gd[1]
    
    # allocate our arrays for the grid quantities
    d.d["ρD"]  = ones(gp[1].dim) # mass density 
    d.d["Φ"]   = ones(gp[1].dim) # gravitational potential

    # allocate particles
    p["x"] = zeros(rank, Npart)
    p["v"] = zeros(rank, Npart)
    p["m"] = 1.0 * length(d.d["ρD"])/Npart # unity for particle mass.

    # define cosmology if it is specified in parameter file
    if haskey(conf, "Cosmology")
        eval(parse(string("global const cosmo=",conf["Cosmology"])))
        if haskey(conf, "InitialRedshift") # set start times from redshift
            zstart = conf["InitialRedshift"]
            u = get_units(cosmo, zstart, zinit=zstart)
            tstart = age_gyr(cosmo, zstart)*(gyr_in_s/u.T)
            zend = conf["FinalRedshift"]
            tend = age_gyr(cosmo, zend)*(gyr_in_s/u.T)
            conf["CurrentTime"] = conf["StartTime"] = tstart
            conf["StopTime"] = tend
            @info("Start Redshift:", zstart, " -> t_start = ",tstart)

            conf["OmegaCDM"] = cosmo.Ω_m
            if haskey(conf,"OmegaBaryon")
                conf["OmegaCDM"] = cosmo.Ω_m - conf["OmegaBaryon"]
            end
            p["m"] .*= conf["OmegaCDM"] # particles
             
        else
            @warn("Specified Cosmology but gave no InitialRedshift! ")
            @warn("Proceed with caution ...")
        end
    end
    
    # call Initializer given as ProblemDefinition in the .conf file
    eval(parse(string("initialvalues=",conf["ProblemDefinition"])))
    
    initialvalues(gp, gd, p)    

    nothing
end

function UniformParticles(gp, gd, p)

    initialize_particles_uniform(p["x"])
    
    nothing
end

function fftfreq(iijk, dims)
    ijk = iijk
    s = div(dims,2)
    k = 2pi .* (ijk-1 - ((ijk .> (s+1)).*dims)) ./dims
end


function PowerSpectrumParticles(gp, gd, p)

    # Initialize a PowerSpectrum if requested
    if haskey(conf, "InputPowerSpectrum")
        eval(parse(string("const ps = ", conf["InputPowerSpectrum"])))
    else
        @critical("Requested to use PowerSpectrumParticles but did not specify InputPowerSpectrum")
    end

    x = p["x"]
    initialize_particles_uniform(x)

    rank = size(x,1)
    dims = conf["ParticleDimensions"]
    ndims = copy(dims)
    ndims[1] = div(dims[1],2)+1
    c = zeros(Complex{Float64},ndims...) #
    Si = zeros(Float64, (rank, (dims[dims .> 1])...))
    for dd in 1:rank
        for k in 1:dims[3], j in 1:dims[2], i in 1:ndims[1]
            kf = fftfreq([i,j,k], dims)
            k2 = dot(kf,kf)
            ka = sqrt(k2)
            gauss = randn(2)
            PSv = sqrt(P(ka, ps))
            ak = Float64(PSv * gauss[1]/k2)
            bk = Float64(PSv * gauss[2]/k2)
#            println(k2,":", PSv)
            if k2 > 0 
                c[i,j,k] = (ak - im * bk)/2 * kf[dd]
            end
        end
        println(summary(c))
        Si[dd,:] = irfft(squeeze(c), dims[1])
#        @show size(Si)
    end

    # Apply displacement field
    for i in 1:size(x,2)5
        for dd in 1:rank
            x[dd,i] += Si[dd,i]
        end
    end
    # Setup velocities : Still need to implement ...
    
    nothing
end


function initialize_particles_uniform(x)
    rank = size(x,1)
    c=1
    if rank==1 
        for i in 1:ParticleDimensions[1]
            x[i] = i
        end
    end
    if rank == 2
        for j in 1:ParticleDimensions[2]
            for i in 1:ParticleDimensions[1]
                x[1,c] = i
                x[2,c] = j
                c += 1
            end
        end
    end
    if rank == 3
        for k in 1:ParticleDimensions[3]
            for j in 1:ParticleDimensions[2]
                for i in 1:ParticleDimensions[1]
                    x[1,c] = i
                    x[2,c] = j
                    x[3,c] = k
                    c += 1 
                end
            end
        end
    end
    x[:] = (x .- 0.5) ./ ParticleDimensions[1:rank]
    nothing
end

function deposit(rho,x,m; interpolation="none")
    for i in eachindex(rho)
        rho[i] = 0.
    end
    
    if interpolation == "none"
        ngpdensity(rho,x,m)
    elseif interpolation == "cic"
        cicdensity(rho,x,m)
    elseif interpolation == "sic"
        sicdensity(rho,x,m)
    end

    
end

function sic_interp_1d(pa, a, x)
    for n in 1:lenth(x)-1
        i = floor(Int64,x[n]*Ngl)+1
        if (x[n+1] > x[n])
            delt = 1
        else
            delt = -1
        end
        cx = x[n]
        inext = i + delt
        xnext = inext/Ngl
        cdx = 
        dxe = abs(x[n+1]-x[n])
        ip1 = (i == Ngl) ? 1 : i+1
        ip1 = (i == 1) ? 1 : i+1
        dx = x[n]*Ngl-i+1
        @inbounds pa[n] = a[1,i]*(1.-dx) + a[1,ip1]*dx
    end
    nothing
end

function ngpdensity(rho,x,m)
    Nglx = convert(Float64,TopgridDimensions[1])
    Ngly = convert(Float64,TopgridDimensions[2])
    Nglz = convert(Float64,TopgridDimensions[3])
    
    if rank == 1
        for i in 1:size(x,2)
            ii = floor(Int64,x[1,i]*Nglx)+1
            rho[ii] += m
        end
    end
    if rank == 2
        for i in 1:size(x,2)
            ii = floor(Int64,(x[1,i]*Nglx)) + 1.
            ij = floor(Int64,(x[2,i]*Ngly)) + 1.
            rho[ii,ij] += m
        end
    end
    if rank == 3
        for i in 1:size(x,2)
            ii = floor(Int64,x[1,i]*Nglx + 1.)
            ij = floor(Int64,x[2,i]*Ngly + 1.)
            ik = floor(Int64,x[3,i]*Nglz + 1.)
            rho[ii,ij,ik] += m
        end
    end
    nothing
end

function cicdensity(rho,x,m)
    Npart = size(x,2)
    Ng = size(rho)
   if rank == 1
        for i in 1:size(x,2)
            ii = floor(Int64,x[1,i]*Ng[1])+1
            ip1 = (ii == Ng[1]) ? 1 : ii+1
            dx = x[1,i]*Ng[1]-ii+1
            rho[ii]  += m*(1-dx)
            rho[ip1] += m*(dx)
        end
    end
    if rank == 2
        for n in 1:Npart
            begin
                ii = floor(Int64,(x[1,n]*Ng[1]))+1
                ip1 = (ii == Ng[1]) ? 1 : ii+1
                dx = x[1,n]*Ng[1]-ii+1
                ij = floor(Int64,(x[2,n]*Ng[2]))+1
                jp1 = (ij == Ng[2]) ? 1 : ij+1
                dy = x[2,n]*Ng[2]-ij+1
                rho[ii,ij]  += m*(1.-dx)*(1.-dy)
                rho[ii,jp1] += m*(1.-dx)*dy
                rho[ip1,ij] += m*dx*(1.-dy)
                rho[ip1,jp1]+= m*dx*dy
            end
        end
    end
    if rank == 3
        for n in 1:Npart
            @inbounds begin
                ii  = floor(Int64,x[1,n]*Ng[1])+1
                ip1 = (ii == Ng[1]) ? 1 : ii+1
                dx  = x[1,n]*Ng[1]-ii+1
                ij  = floor(Int64,x[2,n]*Ng[2])+1
                jp1 = (ij == Ng[2]) ? 1 : ij+1
                dy  = x[2,n]*Ng[2]-ij+1
                ik  = floor(Int64,x[3,n]*Ng[3])+1
                kp1 = (ik == Ng[3]) ? 1 : ik+1
                dz  = x[3,n]*Ng[3]-ik+1
                rho[ii,ij,ik]  += m*(1-dx)*(1-dy)*(1-dz)
                rho[ii,jp1,ik] += m*(1-dx)*dy*(1-dz)
                rho[ip1,ij,ik] += m*dx*(1-dy)*(1-dz)
                rho[ip1,jp1,ik]+= m*dx*dy*(1-dz)
                rho[ii,ij,kp1]  += m*(1-dx)*(1-dy)*dz
                rho[ii,jp1,kp1] += m*(1-dx)*dy*dz
                rho[ip1,ij,kp1] += m*dx*(1-dy)*dz
                rho[ip1,jp1,kp1]+= m*dx*dy*dz
            end
        end
    end
    nothing
end

function make_periodic(x)
    for i in eachindex(x)
       (x[i] < 0.) ? x[i] += 1 : nothing
       (x[i] > 1.) ? x[i] -= 1 : nothing
    end
    nothing
end

function compute_potential(phi,rho,rt)
    if ndims(rho) == 1 potential_1d(phi,rho,rt) end
    if ndims(rho) == 2 potential_2d(phi,rho,rt) end
    if ndims(rho) == 3 potential_3d(phi,rho,rt) end
end

function potential_1d(phi,rho,rt)
    Nglx = TopgridDimensions[1]

    rt[:] = rfft(rho)

    for i in 1:size(rt,1)
        kx = 2pi*(i-1)/Ngl
        Ginv = -4.*(sin(kx/2)^2) # unnormalized
        #            Ginv = -(kx^2+ky^2) # unnormalized
        if (Ginv != 0.0)
            rt[i] /= Ginv
        end
    end
    rt[1] = 0. # make sure singularity is zero
    phi[:] = irfft(rt,size(rho,1)) .* 1/Nglx^2
    nothing
    
end


function potential_2d(phi,rho,rt)
    Nglx = TopgridDimensions[1]
    Ngly = TopgridDimensions[2]

    Ngl2y = div(Ngly,2)

    rt[:] = rfft(rho)

    for j in 1:size(rt,2)
        ky = 2pi*(j-1)/Ngly
        if (j > Ngl2y+1)
            ky = 2pi*(j-Ngly-1)/Ngly
        end
        for i in 1:size(rt,1)
            kx = 2pi*(i-1)/Nglx
            Ginv = -4.*(sin(kx/2)^2+sin(ky/2)^2) # unnormalized
            #            Ginv = -(kx^2+ky^2) # unnormalized
            if !(Ginv == 0.0)
                rt[i,j] /= Ginv
            end
        end
    end
    rt[1,1] = 0. # make sure singularity is zero
    phi[:] = irfft(rt,size(rho,1)) .* 1/Nglx^2
    nothing
    
end



function potential_3d(phi,rho,rt)
    Nglx = TopgridDimensions[1]
    Ngly = TopgridDimensions[2]
    Nglz = TopgridDimensions[3]
    
    Ngl2y = round(Int64, Ngly/2)
    Ngl2z = round(Int64, Nglz/2)

    rt[:]   = rfft(rho)
    for k in 1:size(rt,3)
        kz = 2pi*(k-1)/Nglz
        if (k > Ngl2z+1)
            kz = 2pi*(k-Nglz-1)/Nglz
        end
        for j in 1:size(rt,2)
            ky = 2pi*(j-1)/Ngly
            if (j > Ngl2y+1)
                ky = 2pi*(j-Ngly-1)/Ngly
            end
            for i in 1:size(rt,1)
                kx = 2pi*(i-1)/Nglx
                Ginv = -4.*(sin(kx/2)^2+sin(ky/2)^2+sin(kz/2)^2) # unnormalized   Ginv = -(kx^2+ky^2+kz^2) 
                if !(Ginv == 0.0)
                    rt[i,j,k] /= Ginv
                end
            end
        end
    end
    rt[1,1,1] = 0. # make sure singularity is zero
    phi[:] = irfft(rt,size(rho,1)) .* 1/Nglx^2 # not tested for non cubic grids
    nothing
    
end

function deriv1(xd, x)
    if ndims(x) == 1 deriv1_1d(xd,x) end
    if ndims(x) == 2 deriv1_2d(xd,x) end
    if ndims(x) == 3 deriv1_3d(xd,x) end
    nothing
end

function deriv1_1d(xd,x)
    for i in 1:size(x,1)
        ip1 = (i == Ngl) ? 1 : i+1
        im1 = (i == 1) ? Ngl : i-1
        xd[1,i] = (x[ip1] - x[im1])/2
    end
    fac = convert(Float64,Ngl)
    for i in eachindex(xd)
        xd[i] *= fac
    end
    nothing
end

function deriv1_2d(xd,x)
    Ng1 = size(x,1)
    Ng2 = size(x,2)
    for j in 1:Ng2
        jp1 = (j == Ng2) ? 1 : j+1
        jm1 = (j == 1) ? Ng2 : j-1
        for i in 1:Ng1
            ip1 = (i == Ng1) ? 1 : i+1
            im1 = (i == 1) ? Ng1 : i-1
            xd[1,i,j] = (x[ip1,j] - x[im1,j])/2
            xd[2,i,j] = (x[i,jp1] - x[i,jm1])/2
        end
    end
    fac = convert(Float64,Ngl)
    for i in eachindex(xd)
        xd[i] *= fac
    end
    nothing
end

"""Gradient on a 3D lattice"""
function deriv1_3d(xd,x)
    Ng1 = size(x,1)
    Ng2 = size(x,2)
    Ng3 = size(x,3)
    for k in 1:Ng3
        kp1 = (k == Ng3) ? 1 : k+1
        km1 = (k == 1) ? Ng3 : k-1
        for j in 1:size(x,2)
            jp1 = (j == Ng2) ? 1 : j+1
            jm1 = (j == 1) ? Ng2 : j-1
            for i in 1:size(x,1)
                ip1 = (i == Ng1) ? 1 : i+1
                im1 = (i == 1) ? Ng1 : i-1
                xd[1,i,j,k] = (x[ip1,j,k] - x[im1,j,k])/2
                xd[2,i,j,k] = (x[i,jp1,k] - x[i,jm1,k])/2
                xd[3,i,j,k] = (x[i,j,kp1] - x[i,j,km1])/2
            end
        end
    end
    fac = [convert(Float64,TopgridDimensions[i]) for i in 1:3]
    for i in size(xd,1), j in 1:size(x,2)
        @inbounds xd[i,j] *= fac[i]
    end
    nothing
end


function deriv2(x)
    if ndims(x) == 1 xd = deriv2_1d(x) end
    if ndims(x) == 2 xd = deriv2_2d(x) end
    if ndims(x) == 3 xd = deriv2_3d(x) end
    xd
end

function deriv2_1d(x)
    xd = similar(x)
    for i in 1:size(x,1)
        ip1 = (i == Ngl) ? 1 : i+1
        im1 = (i == 1) ? Ngl : i-1
        xd[i] = x[ip1] + x[im1] - 2*x[i]
    end
    xd*Ngl^2 # 
end

function deriv2_2d(x)
    xd = similar(x)
    for j in 1:size(x,2)
        jp1 = (j == Ngl) ? 1 : j+1
        jm1 = (j == 1) ? Ngl : j-1
        for i in 1:size(x,1)
            ip1 = (i == Ngl) ? 1 : i+1
            im1 = (i == 1) ? Ngl : i-1
            xd[i,j] = x[ip1,j] + x[im1,j] + x[i,jp1] + x[i,jm1] - 4*x[i,j]
        end
    end
    xd*Ngl^2 # 
end

function deriv2_3d(x)
    xd = similar(x)
    for k in 1:size(x,3)
        kp1 = (k == Ngl) ? 1 : k+1
        km1 = (k == 1) ? Ngl : k-1
        for j in 1:size(x,2)
            jp1 = (j == Ngl) ? 1 : j+1
            jm1 = (j == 1) ? Ngl : j-1
            for i in 1:size(x,1)
                ip1 = (i == Ngl) ? 1 : i+1
                im1 = (i == 1) ? Ngl : i-1
                xd[i,j,k] = x[ip1,j,k] + x[im1,j,k] + x[i,jp1,k] + x[i,jm1,k] + x[i,j,kp1] + x[i,j,km1] - 6*x[i,j,k]
            end
        end
    end
    xd*Ngl^2 # divide by h^2
end


function compute_acceleration(a, phi)
    deriv1(a,phi)
    for i in eachindex(a)
        a[i] = -a[i]
    end
    nothing
end

"""For a vector field look up values at location of points)"""
function interpVecToPoints(pa, a, x; interpolation="none")
    if interpolation=="none" # NGP interpolation
        Ng = size(a)[2:end]
        for n in 1:size(x,2)
            i = floor(Int64,x[1,n]*Ng[1])+1
            if rank > 1  j = floor(Int,x[2,n]*Ng[2])+1 end
            if rank > 2  k = floor(Int,x[3,n])*Ng[3]+1 end
            if rank == 1  pa[1,n] = a[1,i]  end
            if rank == 2  pa[:,n] = a[1:2,i,j]  end
            if rank == 3  pa[:,n] = a[1:3,i,j,k]  end
        end
    end
    if interpolation=="cic" # Cloud in Cell interpolation
        if rank == 1
            cic_interp_1d(pa, a, x)
        elseif rank == 2
            cic_interp_2d(pa, a, x)
        elseif rank == 3
            cic_interp_3d(pa, a, x)
        end
    end
    nothing    
end
    
function cic_interp_1d(pa, a, x)
    for n in eachindex(x)
        i = floor(Int64,x[n]*Ngl)+1
        ip1 = (i == Ngl) ? 1 : i+1
        dx = x[n]*Ngl-i+1
        @inbounds pa[n] = a[1,i]*(1.-dx) + a[1,ip1]*dx
    end
    nothing
end


function cic_interp_2d(pa, a, x)
    Ng1 = size(a,2)
    Ng2 = size(a,3)
    for n in 1:size(x,2)
        i = convert(Int64,floor(x[1,n]*Ng1))+1
        ip1 = (i == Ng1) ? 1 : i+1
        dx = x[1,n]*Ng1-i+1
        j = convert(Int64,floor(x[2,n]*Ng2))+1
        jp1 = (j == Ng2) ? 1 : j+1
        dy = x[2,n]*Ng2-j+1
        for m in 1:2
            pa[m,n] = (a[m,i,j]*(1.-dx)*(1.-dy) + a[m,ip1,j]*dx*(1.-dy) +
                       a[m,i,jp1]*(1.-dx)*dy    + a[m,ip1,jp1]*dx*dy)
        end
    end
    nothing
end

function cic_interp_3d(pa, a, x)
    Ng1 = size(a,2)
    Ng2 = size(a,3)
    Ng3 = size(a,4)
    
    for n in 1:size(x,2)
        i = floor(Int64,x[1,n]*Ng1)+1
        ip1 = (i == Ng1) ? 1 : i+1
        dx = x[1,n]*Ng1-i+1
        j = floor(Int64,x[2,n]*Ng2)+1
        jp1 = (j == Ng2) ? 1 : j+1
        dy = x[2,n]*Ng2-j+1
        k = floor(Int64,x[3,n])*Ng3+1
        kp1 = (k == Ng3) ? 1 : k+1
        dz = x[3,n]*Ng3-k+1
        for m in 1:3
            @inbounds pa[m,n] = (a[m,i,j,k]*(1.-dx)*(1.-dy)*(1.-dz) + 
                                 a[m,ip1,j,k]*dx*(1.-dy)*(1.-dz) +
                                 a[m,i,jp1,k]*(1.-dx)*dy*(1.-dz) + 
                                 a[m,ip1,jp1,k]*dx*dy*(1.-dz) +       
                                 a[m,i,j,kp1]*(1.-dx)*(1.-dy)*dz + 
                                 a[m,ip1,j,kp1]*dx*(1.-dy)*dz +    
                                 a[m,i,jp1,kp1]*(1.-dx)*dy*dz +    
            a[m,ip1,jp1,kp1]*dx*dy*dz)
        end
    end
    nothing
end

function drift(x,v,dt)
    for n in eachindex(x)
        x[n] = x[n] + dt*v[n]
    end
    nothing
end


function kick(v,a,dt)
    for n in eachindex(v)
        v[n] = v[n] + dt*a[n]
    end
    nothing
end

#cosmology case
function kick(v,pa,dt,dadt)
#semi implicit integration
    coef = 0.5*dadt*dt
    coef1 = 1.0 - coef
    coef2 = 1.0 / (1.0 + coef)

    for n in eachindex(v)
        v[n] = (coef1*v[n] + dt*pa[n])*coef2
    end
    nothing
end



function evolvePM(gp, gd, p)
    @debug("Entered evolvePM")
    c = conf # convenience

    if isdefined(:cosmo)
        @warn("evolvePM: You defined a cosmology. This Evolve routine will not use it.")
        @warn("evolvePM: You may want to use evolvePMCosmology")
    end

    g = gp[1] # only supporting one grid patch for now
    x = p["x"]
    v = p["v"]
    m = p["m"]
    rho = gd[1].d["ρD"]

    φ = gd[1].d["Φ"]          # potential
    acc   = zeros((3,g.dim...)) 
    
    pa = zeros(x)          # particle accelerations
    rt = rfft(φ) # initialize buffer for real-to complex transform

    dx = 1. /maximum(g.dim) # smallest dx
    dt = InitialDt
    
    c["CurrentTime"] = c["StartTime"]
    c["CurrentCycle"] = c["StartCycle"]
    
    while c["CurrentTime"] < StopTime && c["CurrentCycle"] < StopCycle
        make_periodic(x)
        deposit(rho,x,m,interpolation=ParticleDepositInterpolation)
        rho_mean = mean(rho)
        rho -= rho_mean # the total sum of the density needs to be zero for a periodic signal
        compute_potential(φ, rho, rt)
        compute_acceleration(acc, φ)
        interpVecToPoints(pa, acc, x, interpolation=ParticleBackInterpolation)

        # LeapFrog step. Some codes combine the two drift steps,
        # but then need to worry about x and v being known at different times. 
        drift(x,v,dt/2)
        kick(v,pa,dt)
        drift(x,v,dt/2)

        c["CurrentTime"] += dt

        dt = CourantFactor*dx/maximum(abs(v)+1e-30)
        dt = minimum([dt, MaximumDt])
        
        if ((c["CurrentTime"]+dt) > c["StopTime"] )
            dt = (1. + 1e-15)*(c["StopTime"] - c["CurrentTime"])
            println("final dt = ", dt)
        end
        c["CurrentCycle"] += 1

        analyzePM(gp, gd, p) # check for output
    end

    @info("Evolved to time:", c["CurrentTime"])

end

function powerSpectrum(gp, gd, p)
    x = p["x"]
    make_periodic(x)
    v = p["v"]
    m = p["m"]
    dims = collect(gp[1].dim) # grid dimensions
    box = gp[1].RE - gp[1].LE # box size (may be different from 1.0^3??)
    rho = zeros(gd[1].d["ρD"])
    deposit(rho, x, m, interpolation=ParticleDepositInterpolation)
    rho_mean = mean(rho)
    delta = rho/rho_mean - 1 # density fluctuation

    deltak2 = abs(fft(delta)).^2 * prod(box) / prod(dims)^2
    
    dk = 2pi / maximum(dims) # bin width
    #println("dk", dk)
    lenk = round(Int, maximum(dims)/2) # k array length
    karray = collect(1:lenk) * dk # k array
    power = zeros(lenk) # power spectrum
    counts = zeros(lenk) # array to keep counts
    # we first add up power for each bin, then divide by counts 
    for k in 1:dims[3], j in 1:dims[2], i in 1:dims[1]
        kf = fftfreq([i,j,k], dims)
        ka = sqrt(dot(kf,kf))
        # find the bin where ka falls
        n = round(Int, ka/dk)
        if n >= 1 && n <= lenk
            power[n] += deltak2[i,j,k]
            counts[n] += 1
        end
    end # END looping over all frequencies
    
    power = power./counts

    karray, power
end

function evolvePMCosmology(gp, gd, p)
    @debug("Entered evolvePM")
    if !isdefined(NoName,:cosmo)
        @warn("evolvePMCosmology: No cosmology defined!")
    end

    c = conf # convenience

    g = gp[1] # only supporting one grid patch for now
    x = p["x"]
    v = p["v"]
    m = p["m"]
    rho = gd[1].d["ρD"]
    φ = gd[1].d["Φ"]          # potential
    acc   = zeros((3,g.dim...)) 
    a = 1.0
    ȧ = 0.0
    
    pa = zeros(x)          # particle accelerations
    rt = rfft(φ) # initialize buffer for real-to complex transform

    dx = 1. /maximum(g.dim) # smallest dx
    dt = InitialDt
    zinit = conf["InitialRedshift"]
    
    c["CurrentTime"]  = c["StartTime"]
    c["CurrentCycle"] = c["StartCycle"]

    
    while c["CurrentTime"] < c["StopTime"] && c["CurrentCycle"] < c["StopCycle"]
        make_periodic(x)
        deposit(rho,x,m,interpolation=ParticleDepositInterpolation)

        compute_potential(φ, rho .- mean(rho), rt)
        compute_acceleration(acc, φ)
        interpVecToPoints(pa, acc, x, interpolation=ParticleBackInterpolation)

        a = a_from_t_internal(cosmo, c["CurrentTime"]+dt/2, zinit, zwherea1=zinit)
        ȧ = dadt(cosmo, a, zwherea1=zinit)
        
        # LeapFrog step. Some codes combine the two drift steps,
        # but then need to worry about x and v being known at different times. 
        a = a_from_t_internal(cosmo, c["CurrentTime"] + dt/2/2, zinit, zwherea1=zinit)
        drift(x,v,dt/2/a)
        a = a_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zwherea1=zinit)
        ȧ = dadt_from_t_internal(cosmo, c["CurrentTime"] + dt, zinit, zhwerea1=zinit)
        kick(v,pa,dt/a, ȧ)
        a = a_from_t_internal(cosmo, c["CurrentTime"] + 3.0*dt/2/2, zinit, zwherea1=zinit)
        drift(x,v,dt/2/a)

       
        c["CurrentTime"] += dt
        c["CurrentRedshift"] = z_from_t_internal(cosmo, c["CurrentTime"], zinit)

        dt = CourantFactor*dx/maximum(abs(v)+1e-30)
        dtFrac = DtFractionOfTime*c["CurrentTime"]
        dt = min(dt, dtFrac, MaximumDt)
        
        if ((c["CurrentTime"]+dt) > c["StopTime"] )
            dt = (1. + 1e-6)*(c["StopTime"] - c["CurrentTime"])
            println("final dt = ", dt)
        end
        c["CurrentCycle"] += 1

        analyzePM(gp, gd, p) # check for output
    end

    @info("Evolved to time:", c["CurrentTime"], " Redshift:", c["CurrentRedshift"])

end


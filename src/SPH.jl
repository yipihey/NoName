#module SPH
#export compute_SPH_densities

using NearestNeighbors
import NearestNeighbors.check_input

"""Cubic Spline Kernel popular in Smoothed Particles Hydrodynamics (1,2,3D form)
   call with
   w_ij = W(r,h,ndim=2)
   for twodimensional kernel function
   where r is the (scalar) distance and h the smoothing length
"""
const SPH_default_dim=2

@inline function W{T<:Real}(r::T,h::T;ndim=SPH_default_dim)
    const norm = [4./3., 40/7pi, 8/pi]
    u = r/h
    w = 0.
    if u<=0.5
        w = 1.-6u^2+6u^3
    elseif u<=1
        w = 2(1-u)^3
    end
    w*norm[ndim]/h^ndim
end
# test normalization
# x = collect(1:10000)/10000
# sum(SPH.W(x,1.,ndim=1) *      diff(x)[1])  # should come out to 1 
# sum(SPH.W(x,1.,ndim=2).* x   *diff(x)[1]*2pi)
# sum(SPH.W(x,1.,ndim=3).*x.^2 *(diff(x)[1])*4pi)
W{T<:Real}(r::Vector{T}, h::T;ndim=SPH_default_dim)  = [W(r[i],h,ndim=ndim) for i in eachindex(r)]


""" Derivative of cubic spline kernel
       
       dw ij = dW(r,h, ndim=3)

       ndim can be 1,2 and 3

       dw_ij = dW(r,h)
       uses ndim = SPH_default_dim
       SPH_default_dim should be set as a constant such as
       const SPH_defaul_dim = 2
"""
@inline function dW{T<:Real}(r::T,h::T;ndim=SPH_default_dim)
    const norm = [4./3., 40/7pi, 8/pi]
    u = r/h
    w = 0.
    if u<=0.5
        w = -12u+18u^2
    elseif u<=1
        w = -6(1-u)^2
    end
    w*norm[ndim]/h^(ndim+1)
end

dW{T<:Real}(r::Vector{T}, h::T;ndim=SPH_default_dim) = [dW(r[i],h,ndim=ndim) for i in eachindex(r)]

using Distances
using NearestNeighbors

"""
       Calculate the Euclidian distance assuming all vector components
       values are within the interval [0..1] and assumed periodic.
       The shortest distance to a point or one of its periodic images is
       calculated.
       So
       
       x = ones(3,4) \n
       y = zeros(3,4) \n
       periodic_euclidean(x,y)\n
       will return\n
       4-element Array{Float64,1}:\n
       0.0
       0.0
       0.0
       0.0

       Since [0,0,0] is equal to [1,1,1] in a periodic lattice. 
       """
function periodic_euclidean(a::AbstractArray, b::AbstractArray)
    ld = abs(b .- a)
    res = zeros(size(a,2))
    for j in 1:size(a,2)
        d = 0.
        for i in 1:size(a,1)
            @inbounds c = (ld[i,j] > .5) ? 1-ld[i,j] : ld[i,j]
            d += c*c
        end
        res[j] = sqrt(d)
    end
    res
end
    
type PeriodicEuclidean <: Metric  end

" Find the periodic images we may want to check with the tree"
function periodic_images{T <: AbstractFloat}(p::AbstractArray{T})
    ndim = size(p,1)
    res= zeros(ndim, 2^ndim)
    n = zeros(ndim, 2)
    for i in 1:ndim
        n[i,1] = p[i]
        n[i,2] = (p[i] > 0.5) ? p[i]-1 : p[i]+1
    end
    ni = zeros(p)
    for m in 1:2^ndim
        bity = digits(m-1,2,ndim)
        for i in 1:ndim 
            res[i,m] = n[i,bity[i]+1]
        end
    end
    res
end


function findN_NN{T <: AbstractFloat}(tree::KDTree{T}, points::AbstractArray{T}, k::Int; periodic=false)
    check_input(tree, points)
    @assert tree.reordered==false
    idxs,dists = knn(tree,points, k, true) # returns sorted
    # if the farthest point si far enoguh to cross box edge also check
    # the periodic images within a search radius equal to that. 
    if periodic
        l = 1.*dists[end]
        npt = periodic_images(points)
        nid = idxs
        for i in 2:size(npt,2) # the first (original point) is already done 
            append!(nid, inrange(tree, npt[:,i], l, false))
        end
        dists = periodic_euclidean(tree.data[:,nid], points)
        spi = sortperm(dists)[1:k]
        nid = nid[spi]
        dists = dists[spi]
    end
    # DOES NOT WORK YET> FOR NOW MAKE SURE THE TREE IS NOT REORDERED
#    if tree.reordered
#        println("reorder")
#        @inbounds for j in 1:length(nid)
#            nid[j] = tree.indices[nid[j]]
#        end
#    end

    nid, dists
end


function compute_smoothing_lengths{T <: AbstractFloat}(h::AbstractArray{T},
                                                      tree::KDTree{T},
                                                      Nnghb)
    for i in 1:Npart
        id, dists = findN_NN(t, x[:,i], Nnghb,periodic=true)
        h[i] = dists[Nnghb]
    end
    nothing
end

function compute_SPH_densities_and_h{T <: AbstractFloat}(rho::AbstractArray{T},
                                                         h::AbstractArray{T},
                                                         tree::KDTree{T},
                                                         m, 
                                                         Nnghb)
    @assert size(m,1) == 1
    ndim = size(tree.data,1)
    np = size(tree.data,2)
    for i in 1:np
        id, dists = findN_NN(tree, tree.data[:,i], Nnghb,periodic=true)
        h[i] = dists[Nnghb]
        rho[i] = 0.
        for j in 1:Nnghb
            rho[i] += m*W(dists[j],h[i],ndim=ndim)
        end
    end
    
    nothing
end

function compute_SPH_densities{T <: AbstractFloat}(rho::AbstractArray{T},
                                                   x::AbstractArray{T},
                                                   tree::KDTree{T},
                                                   m, 
                                                   Nnghb)
#    @assert size(m,1) == 1
#    @assert size(tree.data,1) = size(x,1)
    ndim = size(tree.data,1)
    np = size(x,2)
    for i in 1:np
        id, dists = findN_NN(tree, x[:,i], Nnghb,periodic=true)
        h[i] = dists[Nnghb]
        rho[i] = 0.
        for j in 1:Nnghb
            rho[i] += m*W(dists[j],h[i],ndim=ndim)
        end
    end
    
    nothing
end


function compute_SPH_pressures(P,rho,entropy)
    # isothermal equation of state hardcoded for now

    for i in eachindex(rho)
        P[i] = entropy[i]*rho[i]^(GAMMA-1) 
    end
end

function soundspeed(P,rho)
    sqrt(GAMMA*P./rho)
end


function SPH_accelerations_and_update_entropy(a,
                                              P,
                                              entropy,
                                              rho,
                                              h,
                                              tree,
                                              m, Nnghb, dt    )
    @assert size(m,1) == 1
    ndim = size(tree.data,1)
    np = size(tree.data,2)
    dx = zeros(ndim)
    for i in 1:np
        dtEntropy = 0.
        id, dists = findN_NN(tree, tree.data[:,i], Nnghb,periodic=true)
        
        a[:,i] = 0.
        p_over_rho2_i = P[i]/rho[i]^2 
        soundspeed_i = soundspeed(P[i], rho[i])
        for k in 2:Nnghb
            j = id[k]
            r = dists[k]
            dW_i = dW(dists[k],h[i],ndim=ndim)
            dW_j = dW(dists[k],h[j],ndim=ndim)
            vdotr = 0.
            for dim in 1:ndim
                dx[dim] = tree.data[dim,i]-tree.data[dim,j]
                if dx[dim] > 0.5
	            dx[dim] -= 1
                elseif dx[dim] < -0.5
	            dx[dim] += 1
                end

                dvx = v[dim,i] - v[dim,j];
		vdotr += dx[dim] * dvx
            end
            hfc_visc = 0
            if vdotr < 0 # add artificial viscosity of the div < 0
                h_ij = (h[i] + h[j])*0.5
                mu_ij = vdotr* h_ij /( r*r  + 0.01*h_ij*h_ij)
                soundspeed_j = soundspeed(P[j], rho[j])
                c_ij = (soundspeed_i + soundspeed_j)*0.5;
                rho_ij = 0.5 * (rho[i] + rho[j]);
                visc = (- ArtBulkViscConst*c_ij*mu_ij
                        + ArtBulkViscConstB*mu_ij*mu_ij)/rho_ij;
                hfc_visc = visc*0.5*m*(dW_i+dW_j)/r

                dtEntropy += 0.5 * hfc_visc * vdotr
            end

            hfc = hfc_visc + m*(p_over_rho2_i*dW_i+P[j]/rho[j]^2*dW_j)/r

            for dim in 1:ndim
                a[dim,i] -= hfc*dx[dim] 
            end


        end
       entropy[i] += dtEntropy*dt
    end

    nothing
end


function evolveSPH(x,v,P,m,dtin,tfinal; tout=-1)

    pa = zeros(x)          # particle accelerations
    prho = zeros(Npart)    # particle densities
    h = zeros(Npart)       #
    tree = KDTree(x,reorder=false) # construct tree for searching
    compute_SPH_densities_and_h(prho, h, tree, m, Ngb)
    entropy = P./prho.^(GAMMA-1) # initialize entropy from Pressure
    
    dt = dtin
    ttime = 0.

    while ttime < tfinal
        make_periodic(x)
        tree = KDTree(x,reorder=false) # construct tree for searching
        compute_SPH_densities_and_h(prho, h, tree, m, Ngb)
        compute_SPH_pressures(P,prho,entropy)

        cs = maximum(soundspeed(P,prho))
        vm = maximum(v)
        # conservative new timestep
        dt = 0.1*minimum(h)/(cs+vm)
        if ((ttime+dt) > tfinal )
            dt = 1.0001*(tfinal - ttime)
            println("final dt", dt)
        end
        println("dt:", dt)
        
        SPH_accelerations_and_update_entropy(pa, P, entropy,
                                                     prho, h,
                                                     tree, m, Ngb, dt)


        # LeapFrog step. Some codes combine the two drift steps 
        drift(x,v,dt/2)
        kick(v,pa,dt)
        drift(x,v,dt/2)

        ttime += dt
    end


    if debug
        println("Evolved to time:", ttime)
    end
    prho, h, P, pa, entropy
end

function sph_particle_output(fname, x, v, rho, P; overwrite=false)
    if overwrite
        f = h5open(fname,"w")
    else
        f = h5open(fname,"r+")
    end
    g = g_create(f, "particles")
    rank = size(x,1)
    g["x"] = x[1,:]
    g["vx"] = v[1,:]
    if rank > 1
        g["y"] = x[2,:]
        g["vy"] = v[2,:]
    end
    if rank > 2
        g["z"] = x[3,:]
        g["vz"] = v[3,:]
    end
    g["rho"] = rho
    g["P"]   = P
    close(f)
end


#const GAMMA = 5./3  # default adiabatic index
#using PyPlot


#x = rand(2,100000)
#t = KDTree(x,reorder=false)
#points = [0.01,0.5]
#Nnghb  = 32

#Npart = size(x,2)
#h = zeros(Npart)
#rho = zeros(h)
#@time compute_smoothing_lengths(h, t, Nnghb)
#@time compute_SPH_densities(rho, h, t, Nnghb)
#@profile compute_smoothing_lengths(h, t, Nnghb)

#clf()
#plot(x[1,:], x[2,:],   "b.")
#plot(x[1,id], x[2,id], "rx")

#end

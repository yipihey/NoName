# Friends of Friends clustering.  (Tom Abel, 1/2016)
module FOF

using NearestNeighbors
using DataStructures

export groups

"""
Friends of Friends Algorithm: for sets of points x, using linking length l and returning groups with more than minlength members in an array of arrays of indeces inot the point data provided. 
    using Gadfly
    x = rand((2,100)) 
    gps = FOF.groups(x,.1, minlength=3) 
    n = 1
    plot(layer(x=x[1,gps[n]], y=x[2,gps[n]], Geom.point,Theme(default_color=colorant"red")), layer(x=x[1,:], y=x[2,:], Geom.point, Theme(default_color=colorant"blue"))  )

"""
function groups(x, l; minlength=1)
    Npart = size(x)[2]
    tree = KDTree(x) # Build tree with point data provided: KDtree is twice as fast as a balltree
    println("FOF: built tree")
    ds = IntDisjointSets(Npart)
    for i in 1:Npart
        idxs = inrange(tree, x[:,i], l, false) # within search radius
        for j in  eachindex(idxs) # 
            union!(ds,i,idxs[j])
        end
        if (num_groups(ds) == 1) # just in case people use too large a linking length don't waste time
            println("FOF: All points were linked. Exiting." )
            break
        end  # in case everything has been joined already exit
    end
    println("FOF: finished grouping")
    
    idxs = find(ds.ranks) # all non-zero ranks are parent particles in groups
    groupid = [ds.parents[idxs[i]] => i  for i in eachindex(idxs)]
    grouplen = Dict{Int,Int}()
    for i in 1:Npart
        if get(groupid, ds.parents[i], 0) > 0
            grouplen[ds.parents[i]] = get(grouplen, ds.parents[i], 0) + 1
        end
    end

    # now we collect the actual particles in the groups of the length we are interested in
    for (k,v) in grouplen
        if (v < minlength)
            delete!(grouplen, k)
        end
    end

    Ngroups = length(grouplen)
    # and provide them in reverse order with the biggest group first
    sid = sort(collect(grouplen), by = tuple -> last(tuple),rev=true)
    grouplo = [sid[i].first => i for i in 1:length(sid)]
    gps = [[] for i in 1:Ngroups]
    for i in 1:Npart
        if get(grouplen, ds.parents[i], 0) > 0
            push!(gps[grouplo[ds.parents[i]]], i)
        end
    end

    println("FOF: Found ", Ngroups, " with ", minlength, " or more points")
    gps # An array of different  sized arrays containing the particle ids in the groups is returned
end

end

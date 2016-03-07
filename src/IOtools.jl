# io routines for this code
module IOtools

export parse_commandline, write_all, load_all
export write_profiling_results, load_profiling_results

using JLD, HDF5, FileIO
using ArgParse, Logging
@Logging.configure(level=DEBUG)

import Base.squeeze
# allow for easy removal of dimensions of length 1
function squeeze(A::AbstractArray)
    s = size(A)
    dims = tuple(find([dim == 1 for dim in s])...)
    return squeeze(A,dims)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--restart", "-r"
        help = "restart flag. The configuration file is assumed to be a restart file."
        action = :store_true
        "--rconf"
        help="Restart .conf file. File to parse to overriding parameters in restart files. Default: restart.conf"
        default = "restart.conf"
        "configuration_file"
        help = "configuration file name"
        required = false
        arg_type = AbstractString
        default = "particle_mesh.conf"
    end

    return parse_args(s)
end

function readyForOutput(c)

    output = false
    # check for output based on NCycle
    if c["OutputEveryNCycle"] > 0 && mod(c["CurrentCycle"], c["OutputEveryNCycle"]) == 0
        output = true
    end
    #check for Output based on time
    if c["OutputDtDataDump"] != 0 &&
        abs(c["CurrentTime"]-c["LastOutputTime"]) >= c["OutputDtDataDump"]
        output = true
    end

    if output
        c["LastOutputTime"] = c["CurrentTime"]
        c["LastOutputCycle"] = c["CurrentCycle"]
    end
    
    return output
end

function outputGlobalParametersText(fname,conf)
    f = open(fname,"w")
    for (k,v) in conf
        s = string(k, " = ", v,"\n")
        write(f,s)
    end
    close(f)
    @debug("Wrote ", fname)
    nothing    
end

function particle_output(fname, p; overwrite=false)
    if overwrite
        f = h5open(fname,"w")
    else
        f = h5open(fname,"r+")
    end
    g = g_create(f, "particles")
    rank = size(p["x"],1)
    g["x"] = collect(p["x"][1,:])
    g["vx"] = collect(p["v"][1,:])
    if rank > 1
        g["y"] = collect(p["x"][2,:])
        g["vy"] = collect(p["v"][2,:])
    end
    if rank > 2
        g["z"] = collect(p["x"][3,:])
        g["vz"] = collect(p["v"][3,:])
    end
    g["m"] = p["m"]
    close(f)
end

function grid_output(fname, gp, gd; overwrite=false)
    if overwrite
        f = h5open(fname,"w")
    else
        f = h5open(fname,"r+")
    end
    for k in keys(gp)
        cgp = gp[k]
        g = g_create(f, string("grid",cgp.id))
        rank = sum(collect(cgp.dim) .> 1)
        if rank == 1 # for 1D sims store cell centers to aid plotting
            gx = (cgp.LE[1]+(collect(1:cpg.dim[1])-0.5) ./collect(cpg.dim[1]) .*
                  (cgp.RE[1]-cgp.LE[1]))
            g["grid_x"] = gx 
        end
        for i in keys(gd[cgp.id].d)
            g[i] = squeeze(gd[cgp.id].d[i])
            # this attribute helps veusz interpret the range of the image
            attrs(g[i])["vsz_range"] = [gp[1].LE[1], gp[1].LE[2], gp[1].RE[1],gp[1].RE[2]]
        end
    end
    close(f)
end

function write_all(fname, conf, gp, gd, p)
    save(File(format"JLD",fname),
         "conf", conf, "grids", gp, "gridData", gd, "particles", p)
    @info("write_all: Written all data to ", fname)
end

function load_all(fname, conf, gp, gd, p)
    co = load(fname, "conf")
    for (i,v) in co conf[i] = v end
    gpi   = load(fname, "grids")
    for (i,v) in gpi gp[i] = v end
    @show gp
    gdi   = load(fname, "gridData")
    for (i,v) in gdi gd[i] = v end
    pi    = load(fname, "particles")
    for (i,v) in pi p[i] = v end
    
    info("load_all: Loaded all data to ", fname)

    nothing
end



function write_profiling_results(fname)
    r = Profile.retrieve();
    f = open(fname, "w")
    serialize(f, r)
    close(f)
end

function read_profiling_results(fname)
    f = open(fname,"w")
    r = deserialize(f);
    #    ProfileView.view(r[1], lidict=r[2])
end

end # module io

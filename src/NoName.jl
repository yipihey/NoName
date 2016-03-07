module NoName
NONAME_SRC_DIR=dirname(Base.source_path())
println(NONAME_SRC_DIR)

export noname, run_noname

using Logging, AppConf, IOtools
@Logging.configure(level=DEBUG)

include("PrivateUtilities.jl")
include("restart.jl")
include("Hydro.jl") # not a module to access same global space
include("cosmology.jl")
include("ParticleMesh.jl")

println("You want to run this most likely from the shell as:")
println(""" "julia --color=yes -i -e "using NoName; gp, gd, p = run_noname();" particle_mesh.conf" """)

function noname()
    
    @debug("in noname()")
    gp = Dict()    # we keep a hash table to keep track of grids
    gd = Dict()    # empty Dict to store Grid Data
    p  = Dict()    # empty Dict to store particle data
    
    # conf["Initializer"] 
    initializeRoutine(gp,gd,p)
    @info("Finished Initializing.")

    # conf["Evolve"] 
    evolveRoutine(gp,gd,p)
    @info("Finished Evolve.")

    write_all(RestartFileName, conf, gp, gd, p)
    @debug("Wrote restart file:", RestartFileName)

    write_run_summary()
    gp, gd, p
end

function write_run_summary()
    if isinteractive()
        println("You can now access all the data in the variables gp, gd, and p.")
        println("All RunTime parameters are in the dictionary NoName.conf")
        println("To run the code again:")
        println("gp, gd, p = run_noname() ;")
    end
end


function run_noname()
    global parsed_args = parse_commandline() # interpret commandline arguments
    @info("Parsed command line args:")
    for (arg,val) in parsed_args
        @info("  $arg  =>  $val")
    end

    parseconf(string(NONAME_SRC_DIR,"/default.conf")) # load default parameters 
    @info("Read ", string(NONAME_SRC_DIR,"/default.conf"))
    global conf = copy(AppConf.conf)

    # allow user customization
    user_conf_file = string(ENV["HOME"],"/.noname/default.conf")
    if isfile(user_conf_file)
        parseconf(user_conf_file)
        @info("Read ", user_conf_file)
        merge!(conf,AppConf.conf)
    end

    parseconf(parsed_args["configuration_file"]) # parse user config file
    @info("Parsed config files.")
    merge!(conf,AppConf.conf)

    # turn config file to variables (not for the faint of heart)
    dictionary_to_variable_names(:conf)

    # restart? Then make sure to use restart as initializer routine
    if parsed_args["restart"]
        conf["Initializer"] = "restart"
    end
        
    eval(parse(string("initializeRoutine=",conf["Initializer"])))
    eval(parse(string("evolveRoutine=",conf["Evolve"])))
    
    done = false
    if haskey(conf,"ProfilingOn")
        if conf["ProfilingOn"]
            Profile.init(delay=0.01) # Initialize Profiler
            @profile gp, gd, p = noname()
            write_profiling_results(ProfilingResultsFile) # currently broken ...
            @debug("Wrote Profiling Results file:", ProfilingResultsFile)    
            done = true
        end
    end

    if !done
        gp, gd, p = noname()
    end
    gp, gd, p
end


end

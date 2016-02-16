function restart(gp, gd, p)
    # file specified on command line is restart file
    fname = parsed_args["configuration_file"] 
    load_all(fname, gp, gd, p)
    if isfile(parsed_args["rconf"]) # load seperate .conf file if requested
        parseconf(parsed_args["rconf"])
    end
    nothing
end

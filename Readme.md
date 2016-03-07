# NoName
## A code for Computational Cosmology Algorithm Exploration
(2/2016-)



# Introduction
NoName is the beginnings of a Computational Cosmology code of so far 
primarily pedagogical value. 
[Tom Abel](http://tomabel.org/) has been writing this as part of graduate
course at Stanford. It is a flexible enough framework that numerous algorithms 
may be implemented. Some choices made here are influenced by the
[enzo](http://enzo-project.org) and [Gadget](http://wwwmpa.mpa-garching.mpg.de/gadget/) 
codes.


## Current Features
- Particle Mesh dynamics of collisionless particles
- Comoving Coordinates
- Restart files


## Getting Started

The code is intended to be run in the run directory as in:
`julia --color=yes -i -e "using NoName; gp, gd, p = run_noname();" particle_mesh.conf`
Here the `-i` makes julia enter interactive mode after the run is finished and
the code prints some instructions which variables hold the results. 

This assumes that the `LOAD_PATH` contains the directory in which to find
`NoName.jl`. 
You can make sure it is always found by adding 
`push!(LOAD_PATH, "/Users/tabel/Research/codes/noname/src")`
(make sure to change this to your path!) in your
`$HOME/.juliarc.jl` 
file. 

## Dependencies
 Add the following packages with `Pkg.add("Logging")` etc...
 Currently for AppConf we should use `Pkg.checkout("AppConf")` to get the latest version that has not been tagged as a release yet.
 
- [Logging](https://github.com/kmsquire/Logging.jl)      for basic and colorful logging
- [AppConf](https://github.com/tmlbl/AppConf.jl)      for configuration file management
- [HDF5](https://github.com/JuliaLang/HDF5.jl)         data output
- [JLD](https://github.com/JuliaLang/JLD.jl)          restart functionality
- [ArgParse](https://github.com/carlobaldassi/ArgParse.jl)     to parse commandline
- [Cosmology](https://github.com/JuliaAstro/Cosmology.jl)    for cosmological comoving coordinates
- [Roots](https://github.com/JuliaLang/Roots.jl)  for basic root finding
- [NearestNeighbors](https://github.com/KristofferC/NearestNeighbors.jl)  for SPH and Friend of friends group finding


## Caveats
+ No kind of parallelism so far
+ Hydro not fully implemented
+ My first larger julia code framework so very likely full of idiosyncracies 

## Configuration files

Basic code behaviour is influenced by the ___.conf files.
There is one in the `src/default.conf` directory which sets up defaults. 
The user may chose to create directory in their home directory `$HOME/.noname/` 
and store a `default.conf` there to match their preferences. This
user defined file will be read just before the one specified on the command line at
invocation `julia ../../src/noname.jl particle_mesh.conf`. 
The one on the command line has the last word. 

## Submodules and functionality
- [Input PowerSpectra](./doc/InputPowerSpectra.md)
- [Cosmology](./doc/Cosmology.md)

## Parameters

```julia
# This file specifies all the default parameters for our code
# Check the Source or some of parameter.info outputfiles for a full list
HydroMethod = MUSCL
dims        = (30, 1, 1) # main grid dimensions (active zones)
Nghosts     = ( 3, 3, 3) # ghost zones on each side in each dimension
γ           = 1.66667    # equation of state adiabatic index
DualEnergyFormalism=true # 
DomainLeftEdge  = [0. ,0., 0.]  # Size of Domain: Only tested 0..1 so far
DomainRightEdge = [1., 1., 1.]  # 
CourantFactor = 0.3

cosm_Ωm = .3 # Omega Matter
cosm_h = 0.7 # Hubble constant in units of 100 km/s/Mpc

CurrentOutputNumber=0             # Count Outputs 
OutputDirectory = "/tmp/noname"   # Where to store the output
OutputSeparateDirectories=true    # Create directories for each output
LastOutputCycle = -1
LastOutputTime  = 0.
OutputEveryNCycle = 0
OutputDtDataDump = 0.

verbose = false
RestartFileName = /tmp/restart.jld  #
ProfilingResultsFile = profiling_results.bin 
InitialDt = 0.0000001
MaximumDt = 1e30
StartTime = 0
CurrentTime = 0
StopTime = 0
StartCycle = 0
CurrentCycle = 0
StopCycle = 1000000

ParticleDimensions = [1, 1, 1] # number of particles along each dimension
ParticleDepositInterpolation = cic # none and cic implemented so far
ParticleBackInterpolation    = cic # none and cic implemented so far
```

using HydroEOS, GridTypes

function EvolveHydroEuler(gp, gd, p)

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

    nothing
end

function SolveHydroEquations!(gp, gd, p, dt)
    
nothing
end
function ComputeHydroTimeStep(gp, gd, p)
    
    1
end

function InitializeHydroSimulation(gp, gd, p)
    
    global const rank = sum([(i > 1) for i in TopgridDimensions]) # infer dimensionality 
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

    dictionary_to_variable_names(d.d) # shortcuts for the variable names

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
    global const eos = EOS(γ)

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

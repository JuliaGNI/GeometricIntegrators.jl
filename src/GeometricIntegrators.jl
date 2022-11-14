module GeometricIntegrators

    using Reexport

    @reexport using GeometricBase
    @reexport using GeometricBase.Config
    @reexport using GeometricBase.Utils
    @reexport using GeometricEquations

    
    include("Solutions.jl")
    @reexport using .Solutions
    include("Discontinuities.jl")
    @reexport using .Discontinuities
    include("Integrators.jl")
    @reexport using .Integrators
    include("Simulations.jl")
    @reexport using .Simulations
    include("Tableaus.jl")
    @reexport using .Tableaus


    include("methods/methods.jl")
    include("methods/rungekutta.jl")
    
    include("methods/list.jl")

    for m in nameof.(methods)
        @eval export $m
    end


    function __init__()
        add_config(:verbosity, 1)
    end

end

__precompile__()

module GeometricIntegrators

    include("Utils.jl")
    using .Utils

    export Config, Utils, Equations, Solutions, Interpolation, Solvers,
           BasisFunctions, NumericalFluxes, Tableaus, Integrators, Simulations,
           Diagnostics

    include("Config.jl")
    @reexport using .Config
    include("CommonFunctions.jl")
    @reexport using .CommonFunctions
    include("Equations.jl")
    @reexport using .Equations
    include("Solutions.jl")
    @reexport using .Solutions
    include("Interpolation.jl")
    @reexport using .Interpolation
    include("Solvers.jl")
    @reexport using .Solvers
    include("BasisFunctions.jl")
    @reexport using .BasisFunctions
    include("Quadratures.jl")
    @reexport using .Quadratures
    include("NumericalFluxes.jl")
    @reexport using .NumericalFluxes
    include("Tableaus.jl")
    @reexport using .Tableaus
    include("Integrators.jl")
    @reexport using .Integrators
    include("Simulations.jl")
    @reexport using .Simulations
    include("Diagnostics.jl")
    @reexport using .Diagnostics


    function __init__()
        add_config(:verbosity, 1)
    end

end

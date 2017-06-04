__precompile__()

module GeometricIntegrators

    include("Utils.jl")
    using .Utils

    export Utils, Equations, Solutions, Interpolation, Solvers, BasisFunctions,
           Tableaus, Integrators, Diagnostics

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
    include("Tableaus.jl")
    @reexport using .Tableaus
    include("Integrators.jl")
    @reexport using .Integrators
    include("Diagnostics.jl")
    @reexport using .Diagnostics

    export Problems

    include("Problems.jl")

end

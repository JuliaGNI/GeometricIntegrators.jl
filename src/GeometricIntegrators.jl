__precompile__()

module GeometricIntegrators

    include("utils/macro_utils.jl")

    export Utils, Solvers, Equations, Integrators, Tableaus

    include("Fields.jl")
    include("Utils.jl")
    @reexport using .Utils
    include("Equations.jl")
    @reexport using .Equations
    include("Solutions.jl")
    @reexport using .Solutions
    include("Interpolation.jl")
    @reexport using .Interpolation
    include("Solvers.jl")
    @reexport using .Solvers
    include("Polynomials.jl")
    @reexport using .Polynomials
    include("Tableaus.jl")
    @reexport using .Tableaus
    include("Integrators.jl")
    @reexport using .Integrators
    include("Problems.jl")
    @reexport using .Problems

end

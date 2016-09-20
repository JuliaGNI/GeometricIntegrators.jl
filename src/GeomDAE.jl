__precompile__()

module GeomDAE

    include("utils/macro_utils.jl")

    export Utils, Solvers, Equations, Integrators, Tableaus

    include("Utils.jl")
    @reexport using .Utils
    include("Equations.jl")
    @reexport using .Equations
    include("Solvers.jl")
    @reexport using .Solvers
    include("Integrators.jl")
    @reexport using .Integrators
    include("Tableaus.jl")
    @reexport using .Tableaus

end

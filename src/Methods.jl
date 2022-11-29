module Methods

    using GeometricBase.Utils
    using GeometricEquations
    using RungeKutta
    using RungeKutta.Tableaus
    using RungeKutta.PartitionedTableaus

    using ..Integrators
    using ..Integrators.VPRK


    include("methods/methods.jl")
    include("methods/projection.jl")
    include("methods/rungekutta.jl")
    include("methods/splitting.jl")
    include("methods/vprk.jl")
    include("methods/vprk_degenerate.jl")
    include("methods/vprk_projected.jl")
    include("methods/dvi.jl")

    include("methods/list.jl")

end

module Extrapolators

    using Documenter: @doc
    using GeometricBase
    using LinearAlgebra

    using GeometricEquations

    import Base: Callable

    export Extrapolation,
           EulerExtrapolation,
           MidpointExtrapolation,
           HermiteExtrapolation
    export extrapolate!

    include("extrapolation/extrapolation.jl")
    include("extrapolation/aitken_neville.jl")
    include("extrapolation/euler.jl")
    include("extrapolation/hermite.jl")
    include("extrapolation/midpoint.jl")

end

module Extrapolators

    using Documenter: @doc
    using GeometricBase
    using LinearAlgebra

    using GeometricEquations

    import Base: Callable
    import GeometricBase: solutionstep!

    export Extrapolation,
           EulerExtrapolation,
           MidpointExtrapolation,
           HermiteExtrapolation
    export NoInitialGuess
    export extrapolate!, solutionstep!

    include("extrapolation/vectorfields.jl")
    include("extrapolation/extrapolation.jl")
    include("extrapolation/aitken_neville.jl")
    include("extrapolation/euler.jl")
    include("extrapolation/hermite.jl")
    include("extrapolation/midpoint.jl")

end

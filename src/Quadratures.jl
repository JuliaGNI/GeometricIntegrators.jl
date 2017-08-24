__precompile__()

module Quadratures

    using ..Utils

    export Quadrature
    export GaussLegendreQuadrature, LobattoLegendreQuadrature,
           GaussChebyshevQuadrature, LobattoChebyshevQuadrature,
           ClenshawCurtisQuadrature
    export weights

    include("quadratures/quadrature.jl")

end

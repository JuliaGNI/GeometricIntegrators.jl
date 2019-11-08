module Quadratures

    using ..Utils

    export Quadrature
    export RiemannQuadratureLeft, RiemannQuadratureRight
    export GaussLegendreQuadrature, LobattoLegendreQuadrature,
           GaussChebyshevQuadrature, LobattoChebyshevQuadrature,
           ClenshawCurtisQuadrature
    export quadrature, weights

    include("quadratures/quadrature.jl")

end

__precompile__()

module Quadratures

    export Quadrature
    export GaussLegendreQuadrature, GaussLobattoQuadrature,
           ClenshawCurtisQuadrature

    include("quadratures/quadrature.jl")

end

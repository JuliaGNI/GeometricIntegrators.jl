__precompile__()

module Quadratures

    export Quadrature
    export GaussLegendreQuadrature, GaussLobattoQuadrature,
           ClenshawCurtisQuadrature
    export weights

    include("quadratures/quadrature.jl")

end

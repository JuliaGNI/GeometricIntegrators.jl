__precompile__()

module Quadratures

    export Quadrature
    export GaussLegendreQuadrature, GaussLobattoQuadrature,
           ClenshawCurtisQuadrature

    include("Quadrature/quadrature.jl")

end

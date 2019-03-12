__precompile__()

module NumericalFluxes

    using ..Quadratures

    export Flux, FluxPath, FluxPathLinear, FluxPathTrigonometric

    export evaluate_l, evaluate_r, derivative_l, derivative_r

    include("fluxes/flux_path.jl")
    include("fluxes/flux_path_linear.jl")
    include("fluxes/flux_path_trigonometric.jl")

    include("fluxes/numerical_flux.jl")

end

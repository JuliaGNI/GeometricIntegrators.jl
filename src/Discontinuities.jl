module Discontinuities

    using QuadratureRules

    export Discontinuity, PathIntegral, PathIntegralLinear, PathIntegralTrigonometric

    export evaluate_l, evaluate_r, derivative_l, derivative_r

    include("discontinuities/path_integral.jl")
    include("discontinuities/path_integral_linear.jl")
    include("discontinuities/path_integral_trigonometric.jl")

    include("discontinuities/discontinuity.jl")

end

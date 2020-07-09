module Interpolation

    using ..Common
    using ..Utils

    export Interpolator, HermiteInterpolation

    include("interpolation/interpolation.jl")
    include("interpolation/hermite_interpolation.jl")

end

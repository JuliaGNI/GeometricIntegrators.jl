
using GeometricIntegrators
using GeometricIntegrators.Solutions
using SafeTestsets
using Test

include("dataseries_tests.jl")
include("timeseries_tests.jl")

@safetestset "Deterministic Solutions                                                         " begin include("deterministic_solutions_tests.jl") end
@safetestset "Stochastic Solutions                                                            " begin include("stochastic_solutions_tests.jl") end

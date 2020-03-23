
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.Stochastic
using GeometricIntegrators.Integrators.SPARK
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Tableaus
using GeometricIntegrators.Utils
using Test

include("coefficients_tests.jl")
include("rk_tableaus_tests.jl")
include("vprk_tableaus_tests.jl")
include("splitting_tableaus_tests.jl")
include("stochastic_tableaus_tests.jl")

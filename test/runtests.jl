
using GeometricIntegrators
using GeometricIntegrators.Utils
using Test

include("test_functions.jl")
include("test_utils.jl")

include("solvers/linear_solvers_tests.jl")
include("solvers/jacobian_tests.jl")
include("solvers/nonlinear_solvers_tests.jl")
include("basis_functions/vandermonde_tests.jl")
include("basis_functions/lagrange_basis_tests.jl")
include("basis_functions/lagrange_polynomials_tests.jl")
include("equations/equations_tests.jl")
include("oscillator.jl")
include("kubo_oscillator.jl")
include("solutions/dataseries_tests.jl")
include("solutions/timeseries_tests.jl")
include("solutions/solutions_tests.jl")
include("tableaus/coefficients_tests.jl")
include("tableaus/tableaus_tests.jl")
include("integrators/integrators_test.jl")
include("integrators/vprk_integrators_test.jl")
include("integrators/spark_integrators_test.jl")
include("integrators/stochastic_integrators_test.jl")

#include("utils/hdf5_tests.jl")

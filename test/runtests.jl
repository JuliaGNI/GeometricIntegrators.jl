
using GeomDAE
using Base.Test

function fx(t, x, fx)
    fx .= x
end

function fq(t, q, p, fq)
    fq .= q
end

function fp(t, q, p, fp)
    fp .= q.^2
end

include("equations/equations_tests.jl")
include("integrators/timeseries_tests.jl")
include("integrators/solutions_tests.jl")
include("integrators/integrators_test.jl")
include("integrators/tableaus_tests.jl")
include("solvers/linear_solvers_tests.jl")
include("solvers/nonlinear_solvers_tests.jl")
#include("utils/hdf5_tests.jl")

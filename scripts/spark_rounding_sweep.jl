# Outcome rates of three SPARK cases whose result depends on rounding.
#
# Each case runs from initial conditions that differ from q₀ = [1, 1] by k = 0, …, K ULPs in one
# component. A case whose outcome changes with k depends on rounding, so a test of it asserts a
# rate across such a sweep, not the outcome at one start. The tests in
# `test/spark/spark_integrators_tests.jl` and `test/verification/spark_convergence_tests.jl` take
# their bounds from the rates this script prints.
#
#   julia --startup-file=no --project=. scripts/spark_rounding_sweep.jl [K]
#
# The script activates the test environment through TestEnv, which must be on the load path.

using TestEnv
TestEnv.activate()

using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d
using LinearAlgebra
using Logging
using Printf

include(joinpath(@__DIR__, "..", "test", "helpers", "verification_utilities.jl"))
include(joinpath(@__DIR__, "..", "test", "helpers", "rounding_sweep.jl"))

const K = isempty(ARGS) ? 24 : parse(Int, ARGS[1])
const params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)

# The error of a 10-step run with Δt = 0.01 against Gauss(8), or :singular.
function run_error(method, q₀; kw...)
    tspan = (0.0, 0.1)
    ode = odeproblem(q₀; timespan = tspan, timestep = 0.01, parameters = params)
    idae = idaeproblem(q₀; timespan = tspan, timestep = 0.01, parameters = params)
    ref = integrate(ode, Gauss(8))
    try
        sol = with_logger(NullLogger()) do
            integrate(idae, method; kw...)
        end
        return relative_maximum_error(sol.q, ref.q)
    catch e
        e isa SingularException && return :singular
        rethrow()
    end
end

# The convergence order over Δt = 1/10, 1/20, 1/40 on [0, 1], or :throws.
function run_order(method, q₀)
    build(Δt) = idaeproblem(q₀; timespan = (0.0, 1.0), timestep = Δt, parameters = params)
    href(prob) = integrate(
        odeproblem(q₀; timespan = (0.0, 1.0), timestep = timestep(prob),
            parameters = params),
        Gauss(8))
    try
        r = with_logger(NullLogger()) do
            estimate_convergence_order(build, method, 1.0 ./ (10 .* 2 .^ (0:2));
                reference = href, errormetric = (s, r) -> relative_maximum_error(s.q, r.q),
                plateau = 5e-13, integrate_options = (verbosity = 0, warn_iterations = 0))
        end
        return r.order
    catch
        return :throws
    end
end

cases = [
    ("VSPARK(SPARKLobattoIIIBIIIA(2)), error",
        q -> run_error(VSPARK(SPARKLobattoIIIBIIIA(2)), q),
        x -> x === :singular ? "singular" : x < 1e-6 ? "< 1E-6" : "≥ 1E-6"),
    ("VSPARK(SPARKLobABD(3)), error",
        q -> run_error(VSPARK(SPARKLobABD(3)), q; verbosity = 0, warn_iterations = 0),
        x -> x === :singular ? "singular" : x < 8e-6 ? "< 8E-6" : "≥ 8E-6"),
    ("SPARKLobattoIIIBIIIA(2), order", q -> run_order(SPARKLobattoIIIBIIIA(2), q),
        x -> x === :throws ? "throws" :
             isnan(x) ? "NaN" :
             abs(x - 2) ≤ 0.4 ? "order ≈ 2" : "other order")
]

for (label, f, classify) in cases
    outcomes = [f(q) for q in perturbed_starts([1.0, 1.0], K)]
    tally = Dict{String, Int}()
    for x in outcomes
        tally[classify(x)] = get(tally, classify(x), 0) + 1
    end
    nums = filter(x -> x isa Float64 && isfinite(x), outcomes)
    println(label, ", ", length(outcomes), " starts")
    for (key, n) in sort!(collect(tally); by = first)
        @printf("  %-12s %3d  (%.0f %%)\n", key, n, 100n / length(outcomes))
    end
    isempty(nums) ||
        @printf("  range of the values: %.3g … %.3g\n", minimum(nums), maximum(nums))
end
println("VERSION = ", VERSION, ", Sys.ARCH = ", Sys.ARCH, ", check_bounds = ",
    Base.JLOptions().check_bounds)

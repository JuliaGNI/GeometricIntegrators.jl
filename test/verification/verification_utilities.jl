# Shared utilities for the numerical-verification test suite.
#
# These are plain functions (no `@testset`), meant to be `include`d inside each
# per-family `@safetestset` module. They empirically verify that an integrator
# achieves its documented convergence order and, for geometric integrators, that
# it preserves energy/invariants over long integration times.

using GeometricIntegrators
using Test

"""
    loglog_slope(dts, errs)

Least-squares slope of `log(errs)` versus `log(dts)`, i.e. the empirical
convergence order of an error sequence obtained at the timesteps `dts`.
"""
function loglog_slope(dts::AbstractVector, errs::AbstractVector)
    @assert length(dts) == length(errs)
    x = log.(dts)
    y = log.(errs)
    n = length(x)
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

"""
    linear_slope(xs, ys)

Ordinary (non-log) least-squares slope of `ys` versus `xs`. Used as a drift
diagnostic for the per-interval energy error.
"""
function linear_slope(xs::AbstractVector, ys::AbstractVector)
    @assert length(xs) == length(ys)
    n = length(xs)
    x̄ = sum(xs) / n
    ȳ = sum(ys) / n
    d = sum((xs .- x̄) .^ 2)
    return d == 0 ? zero(ȳ) : sum((xs .- x̄) .* (ys .- ȳ)) / d
end

# Default error metric: relative maximum error in the position component `q`.
default_errormetric(sol, ref) = relative_maximum_error(sol, ref).q

"""
    estimate_convergence_order(problem_builder, method, timesteps;
                               reference, errormetric, plateau)

Integrate `method` at each `Δt` in `timesteps` (building the problem via
`problem_builder(Δt)`), measure the error against `reference(prob)` with
`errormetric(sol, ref)`, and return a NamedTuple `(order, dts, errs, mask)` where
`order` is the log-log slope over the points that lie above the roundoff plateau.

* `problem_builder :: Δt -> problem`
* `reference       :: problem -> solution`   (e.g. `exact_solution` or a high-order run)
* `errormetric     :: (sol, ref) -> Real`
* `plateau`        : points with `err ≤ plateau`, or non-decreasing errors, are discarded
"""
function estimate_convergence_order(problem_builder, method, timesteps;
        reference,
        errormetric = default_errormetric,
        plateau = 1e-13)
    dts = collect(float.(timesteps))
    errs = map(dts) do Δt
        prob = problem_builder(Δt)
        errormetric(integrate(prob, method), reference(prob))
    end
    # Keep points above the roundoff floor that are still strictly decreasing;
    # this drops the flattened tail where high-order methods hit machine precision.
    mask = trues(length(errs))
    for i in eachindex(errs)
        mask[i] = errs[i] > plateau && (i == 1 || errs[i] < errs[i-1])
    end
    ord = count(mask) >= 2 ? loglog_slope(dts[mask], errs[mask]) : NaN
    return (order = ord, dts = dts, errs = errs, mask = mask)
end

"""
    test_convergence_order(problem_builder, method, timesteps; kwargs...)

Assert that the empirical convergence order of `method` matches `expected`
(default `order(method)`) within `atol`, using at least `minpoints` valid points.
Returns the result NamedTuple from [`estimate_convergence_order`](@ref).
"""
function test_convergence_order(problem_builder, method, timesteps;
        reference,
        errormetric = default_errormetric,
        expected = order(method),
        atol = 0.35,
        plateau = 1e-13,
        minpoints = 3,
        label = string(nameof(typeof(method))))
    @assert expected isa Number "pass an explicit `expected` order for $(label)"
    res = estimate_convergence_order(problem_builder, method, timesteps;
        reference, errormetric, plateau)
    @testset "$(label): order ≈ $(expected)" begin
        @test count(res.mask) >= minpoints
        @test isapprox(res.order, expected; atol = atol)
    end
    return res
end

"""
    energy_error(sol, hamiltonian, params; partitioned)

Relative energy error time series `(H(t) - H(0)) / H(0)` computed with the
problem's `hamiltonian`. Uses the partitioned `(t, q, p, params)` signature when
`partitioned` is true, otherwise the `(t, q, params)` signature.
"""
function energy_error(sol, hamiltonian, params; partitioned::Bool)
    _, errds = partitioned ?
        compute_invariant_error(sol.t, sol.q, sol.p, params, hamiltonian) :
        compute_invariant_error(sol.t, sol.q, params, hamiltonian)
    return errds
end

"""
    test_energy_behaviour(prob, method; hamiltonian, partitioned,
                          interval_length, energy_tol, drift_tol, label)

Integrate a long trajectory and assert (a) the relative energy error stays
bounded by `energy_tol`, and (b) the per-interval maximum energy error exhibits
no secular growth (linear drift slope below `drift_tol`). `nt` (from `prob`) must
be divisible by `interval_length`.
"""
function test_energy_behaviour(prob, method; hamiltonian, partitioned::Bool,
        interval_length = 100, energy_tol, drift_tol,
        label = string(nameof(typeof(method))))
    sol = integrate(prob, method)
    params = parameters(prob)
    errds = energy_error(sol, hamiltonian, params; partitioned)
    tdrift, idrift = compute_error_drift(sol.t, errds, interval_length)
    @testset "$(label): energy behaviour" begin
        @test maximum(abs.(errds.d)) < energy_tol
        @test abs(linear_slope(collect(tdrift.d), collect(idrift.d))) < drift_tol
    end
    return (errds = errds, tdrift = tdrift, idrift = idrift)
end

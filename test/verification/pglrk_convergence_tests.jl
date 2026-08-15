using GeometricIntegrators
using GeometricIntegratorsBase
using GeometricProblems.HarmonicOscillator
using GeometricProblems.LotkaVolterra2d
using Logging
using Test

include("verification_utilities.jl")

# Projected Gauss-Legendre Runge-Kutta: the energy-preserving Gauss collocation method
# of Brugnano, Iavernaro & Trigiante (SINUM 50(6), 2012).
#
# Two things are verified, on the *nonlinear* Lotka-Volterra problem where they are
# actually distinguishable:
#
#   1. The order is that of the underlying Gauss tableau (2s). This is expected because
#      the skew perturbation satisfies bᵀA = 0 and A·1 = 0, so a(λ) retains the
#      quadrature and row-sum conditions for every λ.
#   2. The energy is conserved to the tolerance of the λ solve, whereas plain Gauss
#      exhibits an O(h^{2s}) energy error. On the harmonic oscillator both are at
#      machine precision — Gauss already conserves quadratic invariants exactly — so
#      the linear problem cannot distinguish them and is not used for the energy test.
#
# `VPRKpTableau` is covered here too rather than in its own file, because it is built on
# the same `CoefficientsPGLRK` tableau family — with the Dirac constraint substituted for
# the energy condition — and the two are best read side by side.

const T = 2.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

# Counts the statuses PGLRK hands to `check_solver_status`, for the testset of the same name below.
# The override only counts and returns its argument, so it leaves every assertion that runs under
# it exactly as it was. Note how far "under it" reaches: a method is global to the *session*, and
# `runtests.jl` drives its files with `@safetestset`, which is a fresh module in the same process —
# so this counts every PGLRK integration in every file that runs after this one too. Nothing reads
# the counter there, and the reset below is immediately before each integration it is read after
# rather than once here, so that is harmless; it is not file-local, though.
const PGLRK_STATUS_CHECKS = Ref(0)

function GeometricIntegratorsBase.check_solver_status(status, int::GeometricIntegrator{<:PGLRK})
    PGLRK_STATUS_CHECKS[] += 1
    status
end

# Run `method` on `prob` and report both of the things the testset below distinguishes routes by:
# how many statuses reached the hook, and how many steps fell back to plain Gauss. The fallback is
# read off the `@debug` line `solve_λ!` emits on that branch, which is why the integration runs
# under a `TestLogger` — it is the only outward sign of which way through `solve_λ!` a step went,
# and the whole point here is to assert the route rather than to assume it.
function count_checks_and_fallbacks(prob, method)
    PGLRK_STATUS_CHECKS[] = 0
    logger = TestLogger(min_level=Logging.Debug)
    with_logger(() -> integrate(prob, method), logger)
    (checks=PGLRK_STATUS_CHECKS[],
     nfallback=count(r -> occursin("falling back to plain Gauss", r.message), logger.logs))
end

build(Δt) = LotkaVolterra2d.odeproblem(q₀; timespan=(0.0, T), timestep=Δt, parameters=params)
ref(prob) = integrate(LotkaVolterra2d.odeproblem(q₀; timespan=timespan(prob), timestep=timestep(prob), parameters=params), Gauss(8))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))

# relative energy error over the whole trajectory
function energy_drift(prob, method)
    sol = integrate(prob, method)
    ps = parameters(prob)
    h = [LotkaVolterra2d.hamiltonian(sol.t[i], sol.q[i], ps) for i in eachindex(sol.q)]
    return maximum(abs.(h .- h[1])) / abs(h[1])
end

@testset "$(rpad("Projected Gauss-Legendre Runge-Kutta convergence",80))" begin

    @testset "convergence order" begin
        # Measured on this problem and step range: PGLRK 5.69 / 7.74 against
        # Gauss 5.88 / 7.87 — i.e. PGLRK matches the underlying method, and both sit a
        # little under the nominal order here, which is a property of the problem and
        # step range rather than of the projection.
        #
        # `atol` is set from those measurements rather than left at the 0.35 default, which
        # would leave 0.04 of headroom at s = 3 and 0.09 at s = 4 — tight enough that BLAS
        # or dependency drift could flip the assertion.
        test_convergence_order(build, PGLRK(3), steps(2, 3); reference=ref, expected=6, atol=0.5, label="PGLRK(3)")
        test_convergence_order(build, PGLRK(4), steps(2, 3); reference=ref, expected=8, atol=0.5, label="PGLRK(4)")
    end

    @testset "energy conservation vs. plain Gauss" begin
        prob = LotkaVolterra2d.odeproblem(q₀; timespan=(0.0, 10.0), timestep=0.1, parameters=params)

        # Measured: PGLRK(3) 3.6E-15 vs Gauss(3) 7.3E-10 — five orders of magnitude, and
        # at the tolerance of the λ bisection rather than at O(h^6).
        pg3 = energy_drift(prob, PGLRK(3))
        g3 = energy_drift(prob, Gauss(3))
        @test pg3 < 4E-15
        @test pg3 < g3 / 1000

        # At s = 4 the plain method's energy error is already near the λ-solve tolerance,
        # so the improvement is smaller; the projection must still not make it worse.
        pg4 = energy_drift(prob, PGLRK(4))
        g4 = energy_drift(prob, Gauss(4))
        @test pg4 < 4E-13
        @test pg4 < 10 * g4
    end

    @testset "the λ bisection's probes do not reach the solver-status hook" begin
        # PGLRK is the one method here whose stage solves are *probes*: `energy_residual!` runs a
        # full nonlinear solve per trial λ, driven by `SimpleSolvers.bisection`, including at the
        # endpoints ±λmax where a solve that struggles is expected rather than exceptional. Only
        # the accepted λ's solve — the one whose result the step keeps — is handed to
        # `check_solver_status`, which is the hook a caller overrides to reject a step that did not
        # converge. Hand it every probe and such a caller rejects bisection probes instead.
        #
        # Nothing asserted that, and nothing would have: the hook's default returns its argument, so
        # calling it ten times too often is invisible to every other test. Count the calls.
        # The failure this guards against is not "the hook is never called" but "the hook is called
        # once per probe", i.e. once per bisection trial rather than once per step — so the
        # assertion is an exact count and not a bound.
        #
        # `solve_λ!` has three ways through it and a separate `check_solver_status` call site on
        # each, so all three are driven here — and *which* one a run took is asserted rather than
        # assumed, from the fallback count, since a hook count of `nsteps` is what every route is
        # supposed to produce and so cannot itself tell them apart.
        nsteps = 5
        lv = LotkaVolterra2d.odeproblem(q₀; timespan=(0.0, nsteps * 0.1), timestep=0.1, parameters=params)

        # Route 1 — the bisection locates a root. On this problem plain Gauss(3) drifts to
        # |h − h₀| ≈ 6E-11 within these five steps, nowhere near the ftol ≈ 3.6E-15 of the early
        # return, so the first probe never satisfies it and the bisection runs; that it also finds
        # a root, rather than falling back, is what `nfallback == 0` says.
        r = count_checks_and_fallbacks(lv, PGLRK(3))
        @test r.nfallback == 0
        @test r.checks == nsteps

        # Route 2 — the bisection runs and finds no sign change in [±λmax], so the step falls back
        # to plain Gauss. Note that λmax does *not* reach the early return, whose condition is on
        # the λ = 0 residual alone: on this problem this is a fallback on every step, not an early
        # return, and `nfallback == nsteps` is what distinguishes it from route 1.
        r = count_checks_and_fallbacks(lv, PGLRK(3; λmax=1e-300))
        @test r.nfallback == nsteps
        @test r.checks == nsteps

        # Route 3 — the early return: the unperturbed method already conserves the energy to
        # tolerance, so `solve_λ!` returns before the bisection. Gauss preserves the quadratic
        # invariant of a linear system exactly, so the harmonic oscillator is where that branch is
        # reachable (|h − h₀| ≈ 7E-18 against ftol ≈ 1.8E-15). Keeping the λmax of route 2 is what
        # makes this an assertion instead of a hope: had the early return *not* been taken, the
        # bisection would have run in [±1E-300] and fallen back on every step exactly as it does
        # there, so `nfallback == 0` here means the bisection never ran at all.
        ho = HarmonicOscillator.odeproblem(; timespan=(0.0, nsteps * 0.1), timestep=0.1)
        r = count_checks_and_fallbacks(ho, PGLRK(3; λmax=1e-300))
        @test r.nfallback == 0
        @test r.checks == nsteps
    end

    @testset "reduces to Gauss when λ = 0" begin
        # With λmax driven to (almost) zero the method must fall back to plain Gauss.
        prob = build(0.25)
        for s in 3:4
            @test relative_maximum_error(integrate(prob, PGLRK(s; λmax=1e-300)).q,
                                        integrate(prob, Gauss(s)).q) < 4E-16
        end
    end

    @testset "VPRKpTableau on a degenerate Lagrangian" begin
        # `VPRKpTableau` uses the same tableau family, with the Dirac constraint
        # ϑ(qₙ₊₁) − pₙ₊₁ = 0 in place of the energy condition. Its point is to recover on
        # a *degenerate* Lagrangian the order that a plain VPRK loses there.
        # The nonlinear solver occasionally runs to its iteration cap on these harder
        # systems and emits benign warnings. Silencing them needs both kwargs:
        # `verbosity` gates the stagnation and line-search messages, `warn_iterations` the
        # "Solver took N iterations" cap message.
        lbuild(Δt) = LotkaVolterra2d.lodeproblem(q₀; timespan=(0.0, T), timestep=Δt, parameters=params)
        emq(sol, r) = relative_maximum_error(sol.q, r.q)

        # Measured on this problem (T = 2, Δt = 1 … 1/8): plain VPRK(Gauss(s)) suffers the
        # documented order reduction to ≈ s on a degenerate Lagrangian (3.6 / 4.1 / 5.8 for
        # s = 3 / 4 / 5), whereas VPRKpTableau reaches 5.3 / 9.5 / 9.8. The projection
        # roughly doubles the achieved order.
        st = T ./ (2 .* 2 .^ (0:3))
        for s in 3:5
            rp = estimate_convergence_order(lbuild, VPRKpTableau(s), st; reference=ref, errormetric=emq, integrate_options=(verbosity=0, warn_iterations=0,))
            rv = estimate_convergence_order(lbuild, VPRK(Gauss(s)), st; reference=ref, errormetric=emq, integrate_options=(verbosity=0, warn_iterations=0,))
            @test rp.order > rv.order + 1
        end

        # A clean order *number*, however, is not reproducible: the multiplier introduces a
        # Δt-independent error floor (≈1E-11 at s = 4, ≈1E-13 at s = 5) and the error
        # sequence is non-monotone, so the fitted slope swings between ≈5.9 and ≈9.5
        # depending on the step window. Tightening the solver tolerances (f_abstol 1E-14,
        # 1E-16, extra Newton iterations) does not change the floor, so it is a property of
        # the method and not of the solve. `VPRKpTableau` has no published order theorem —
        # see its docstring — so this is recorded rather than asserted.
        @test_broken all(
            isapprox(estimate_convergence_order(lbuild, VPRKpTableau(s), T ./ (4 .* 2 .^ (0:4)); reference=ref, errormetric=emq, integrate_options=(verbosity=0, warn_iterations=0,)).order, 2s; atol=0.35)
            for s in 4:5
        )
    end

    @testset "harmonic oscillator (linear reference)" begin
        # Here the exact solution is available; both PGLRK and Gauss are at machine
        # precision, so this only guards against regressions in the ODE path.
        hobuild(Δt) = HarmonicOscillator.odeproblem(; timespan=(0.0, 1.0), timestep=Δt)
        test_convergence_order(hobuild, PGLRK(3), 1.0 ./ (4 .* 2 .^ (0:2)); reference=exact_solution, expected=6, minpoints=2, label="PGLRK(3) harmonic")
    end

end

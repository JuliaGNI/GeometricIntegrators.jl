using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d
using LinearAlgebra
using RungeKutta
import GeometricIntegratorsBase
using Test

const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt * nt)

ode = odeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
hdae = hdaeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
idae = idaeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
pdae = pdaeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
ldae = ldaeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
ldae_slrk = ldaeproblem_slrk(q₀; timespan=tspan, timestep=Δt, parameters=params)

ref = integrate(ode, Gauss(8))

# Several known-broken / order-reduced SPARK methods below genuinely diverge, hit
# singular stage systems, or reduce order (see docs/src/audit.md, third pass).
# Their solves emit nonlinear-solver stagnation, iteration-cap and backtracking
# line-search warnings that are correct symptoms and not fixable via the solver (a
# different solver, line search or iteration cap does not help). Passing both
# `verbosity = 0` and `warn_iterations = 0` silences those (correct) warnings for the
# affected integrations: under SimpleSolvers 0.10 the stagnation and line-search
# messages are gated on `verbosity`, while the "Solver took N iterations" cap message
# is gated on `warn_iterations`. Since GeometricIntegratorsBase merges these into the
# method's `default_options`, the tolerances and `min_iterations` are unchanged and only
# logging is affected, so the measured errors are unaffected.

# Where the mechanism is unambiguous, the known-broken cases below assert *why* each fails
# rather than merely that it is inaccurate: a structurally singular stage system throws a
# `SingularException` (`@test_throws`); an order-deficient/unstable method whose per-step
# solves are healthy is asserted `converged`; a solve that stalls at the residual floor is
# asserted `stalled`. `nlsolve_outcome` drives one step manually and reads the nonlinear
# solver's status (SimpleSolvers 0.10 exports these, but not by name). The remaining cases
# are only *marginally* singular — they converge or zero-pivot depending on rounding — so
# neither assertion is reliable and they stay `@test_broken`.
using GeometricIntegratorsBase: solver, solverstate
using SimpleSolvers: status, isconverged, isstalled, config
function nlsolve_outcome(prob, method)
    int = GeometricIntegrator(prob, method; verbosity=0, warn_iterations=0)
    ss = GeometricIntegratorsBase.solutionstep(int, GeometricIntegratorsBase.initialstate(prob))
    GeometricIntegratorsBase.reset!(ss, Δt)
    GeometricIntegratorsBase.integrate!(ss, int)
    st = status(solver(int), solverstate(int))
    (converged = isconverged(st), stalled = isstalled(st, config(solver(int))))
end


@testset "$(rpad("SLRK integrators",80))" begin

    sol = integrate(ldae_slrk, SLRKLobattoIIIAB(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIAB(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIAB(4))
    @test relative_maximum_error(sol.q, ref.q) < 8E-16


    sol = integrate(ldae_slrk, SLRKLobattoIIIBA(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIBA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIBA(4))
    @test relative_maximum_error(sol.q, ref.q) < 8E-16


    sol = integrate(ldae_slrk, SLRKLobattoIIICC̄(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIICC̄(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIICC̄(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIIC̄C(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIC̄C(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIC̄C(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIID(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIID(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIID(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-7

    sol = integrate(ldae_slrk, SLRKLobattoIIIE(3))
    @test relative_maximum_error(sol.q, ref.q) < 8E-12

    sol = integrate(ldae_slrk, SLRKLobattoIIIE(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15

end


@testset "$(rpad("SLRK stage-Jacobian conditioning (finding S10)",80))" begin

    # The null-vector multiplier μ must enter the primary-constraint residual row only.
    # Until the fifth pass it was added to the momentum-stage row as well; since that
    # row lives in Z-space (P = p + h·Z), the same coefficient there carries an extra
    # factor h, the two contributions combine to (1-h)·μ·dᵢ/b̄ᵢ, and the μ column of
    # the stage Jacobian drops out of the reduced system at Δt = 1 — making it exactly
    # singular there and ill-conditioned nearby.
    #
    # The integrator tests above run at Δt = 0.01, where this is invisible, and the
    # solver still converges at Δt = 0.99 either way, so only the conditioning of the
    # stage Jacobian discriminates. With this central-difference stencil, measured over
    # the twelve cases below:
    #
    #             |      Δt = 0.99      |        Δt = 1
    #   with bug  |  1.2e3 … 6.5e4      |  7.7e10 … 2.0e12
    #   fixed     |  2.1e1 … 1.6e3      |  2.1e1 … 2.6e3
    #
    # So it is the Δt = 1 row that carries the regression — every one of the twelve
    # breaches the 1e4 bound there — while at Δt = 0.99 only two of them do
    # (SLRKLobattoIIIC̄C(2) at 2.0e4 and SLRKLobattoIIICC̄(3) at 6.5e4). The Δt = 0.99
    # assertions are kept because the 1/(1-Δt) growth is the signature, but they are not
    # what makes the test fail. Worst case for the fixed code anywhere on
    # Δt ∈ {0.1, 0.5, 0.9, 0.99, 1.0} is 2.6e3, i.e. the bound clears it by a factor of
    # four. (The 1e11 is the stencil's noise floor standing in for an exactly singular
    # matrix — step 2 of `scripts/slrk_verification.jl` differentiates exactly and gets
    # 3.1e17. Either way the bound separates fixed from buggy by seven orders of
    # magnitude at Δt = 1, so the cheap stencil is enough for a regression test.)

    "central-difference stage Jacobian of `residual!` at the initial guess"
    function stage_jacobian(method, Δt)
        prob = ldaeproblem_slrk(q₀; timespan=(t₀, t₀ + 10Δt), timestep=Δt, parameters=params)
        int = GeometricIntegrator(prob, method)
        solstep = GeometricIntegratorsBase.solutionstep(int, GeometricIntegratorsBase.initialstate(prob))
        GeometricIntegratorsBase.reset!(solstep, Δt)

        sol = GeometricIntegratorsBase.current(solstep)
        prm = GeometricIntegratorsBase.parameters(solstep)
        GeometricIntegratorsBase.initial_guess!(sol, GeometricIntegratorsBase.history(solstep), prm, int)

        x = copy(GeometricIntegratorsBase.nlsolution(int))
        n = length(x)
        J = zeros(n, n)
        rp = zeros(n)
        rm = zeros(n)
        ε = 1e-7
        for j in 1:n
            xp = copy(x); xp[j] += ε
            xm = copy(x); xm[j] -= ε
            GeometricIntegratorsBase.residual!(rp, xp, sol, prm, int)
            GeometricIntegratorsBase.residual!(rm, xm, sol, prm, int)
            J[:, j] .= (rp .- rm) ./ (2ε)
        end
        J
    end

    for ctor in (SLRKLobattoIIIAB, SLRKLobattoIIIBA, SLRKLobattoIIICC̄,
                 SLRKLobattoIIIC̄C, SLRKLobattoIIID, SLRKLobattoIIIE), s in (2, 3)
        for Δt in (0.99, 1.0)
            @test cond(stage_jacobian(ctor(s), Δt)) < 1E4
        end
    end

end



@testset "$(rpad("SPARK solver size (finding S18)",80))" begin

    # Two invariants that were both broken and both invisible.
    #
    # First: `solversize` must be the *same* generic function the rest of the package
    # extends. `src/SPARK.jl` did not import it from GeometricIntegratorsBase, so the ten
    # definitions in `src/spark/` built a second, unrelated function under the same name —
    # GeometricIntegratorsBase.solversize covered every non-SPARK method and raised a
    # MethodError on every SPARK one.
    #
    # Second: `solversize` must be the full length of the solver's unknown vector. The
    # null-vector multiplier μ contributes D unknowns, and that term used to be added by
    # the cache constructor rather than by `solversize`, so the length was defined in two
    # files at once with nothing checking they agreed. Before the fix the two differed by
    # exactly D for every method carrying a null vector — e.g. 16 against 18 for
    # SLRKLobattoIIID(2) and for TableauVSPARKLobattoIIIAB(2).

    @test GeometricIntegrators.SPARK.solversize === GeometricIntegratorsBase.solversize

    # one method per family; both null-vector branches, and the `true` branch on two
    # different families so it is not just an SLRK property
    solver_size_cases = (
        (idae,      SPARKGLRK(1),                          false),
        (idae,      VSPARK(SPARKGLRKLobattoIIIAIIIB(1)),   false),
        (idae,      TableauVSPARKGLRKpMidpoint(1),          false),
        (idae,      TableauVSPARKGLRKpSymmetric(2),         false),
        (idae,      TableauVSPARKLobattoIIIAB(2),           true),
        (pdae,      TableauHPARKGLRK(1),                    false),
        (ldae_slrk, SLRKLobattoIIID(2),                     true),
        (ldae_slrk, SLRKLobattoIIIAB(3),                    true),
    )

    for (prob, method, hasnull) in solver_size_cases
        int = GeometricIntegrator(prob, method)
        @test GeometricIntegrators.SPARK.hasnullvector(method) == hasnull
        @test GeometricIntegratorsBase.solversize(prob, method) ==
              length(GeometricIntegratorsBase.nlsolution(int))
    end

end



@testset "$(rpad("SPARK integrators",80))" begin

    sol = integrate(idae, SPARKGLRK(1))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, SPARKGLRK(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11


    # TODO: Check Errors !!!  GLVPRK should do much better ! (maybe problem with R∞?)

    sol = integrate(idae, SPARKGLVPRK(1))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, SPARKGLVPRK(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6


    sol = integrate(idae, SPARKLobABC(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, SPARKLobABC(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, SPARKLobABC(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(idae, SPARKLobABD(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, SPARKLobABD(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, SPARKLobABD(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(idae, SPARKLobattoIIIAIIIB(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, SPARKLobattoIIIAIIIB(4))
    @test relative_maximum_error(sol.q, ref.q) < 8E-10


    sol = integrate(idae, SPARKGLRKLobattoIIIAIIIB(1))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(idae, SPARKGLRKLobattoIIIAIIIB(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    # order reduction: accuracy plateaus ~2E-4, does not improve with s
    sol = integrate(idae, SPARKGLRKLobattoIIIAIIIB(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-4


    sol = integrate(idae, SPARKGLRKLobattoIIIBIIIA(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(idae, SPARKGLRKLobattoIIIBIIIA(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    # order reduction: accuracy plateaus ~2E-4, does not improve with s
    sol = integrate(idae, SPARKGLRKLobattoIIIBIIIA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-4


    # --- known-broken (see docs/src/audit.md) ---
    # These are @test_broken rather than mechanism assertions because their stage systems
    # are only *marginally* singular: the LU either converges to an order-deficient answer
    # or zero-pivots depending on rounding, so neither `@test`(converged) nor `@test_throws`
    # is reliable — @test_broken robustly tolerates both.

    # order-deficient at s=2 (only order 1, meas 0.18)
    @test_broken relative_maximum_error(integrate(idae, SPARKLobattoIIIAIIIB(2); verbosity=0, warn_iterations=0).q, ref.q) < 1E-6

    # SPARKLobattoIIIBIIIA diverges (non-finite at s=2; finite-but-wrong at s=3/s=4, meas 0.58 / 0.03)
    @test_broken relative_maximum_error(integrate(idae, SPARKLobattoIIIBIIIA(2); verbosity=0, warn_iterations=0).q, ref.q) < 1E-6
    @test_broken relative_maximum_error(integrate(idae, SPARKLobattoIIIBIIIA(3); verbosity=0, warn_iterations=0).q, ref.q) < 1E-6
    @test_broken relative_maximum_error(integrate(idae, SPARKLobattoIIIBIIIA(4); verbosity=0, warn_iterations=0).q, ref.q) < 2E-10

end


@testset "$(rpad("VPARK integrators",80))" begin

    sol = integrate(idae, TableauSymplecticProjection(:glrk1ps, TableauGauss(1), TableauGauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauSymplecticProjection(:glrk2ps, TableauGauss(2), TableauGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauGausspSymplectic(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauGausspSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauLobattoIIIAIIIBpSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, TableauLobattoIIIAIIIBpSymplectic(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, TableauLobattoIIIAIIIBpSymplectic(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15

    sol = integrate(idae, TableauLobattoIIIBIIIApSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauLobattoIIIBIIIApSymplectic(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-11

    sol = integrate(idae, TableauLobattoIIIBIIIApSymplectic(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


# TODO: Replace idae with ldae !!!

@testset "$(rpad("VSPARK integrators",80))" begin

    # converges to 2.7E-10 but the backtracking line search struggles in the tail
    # (many "did not satisfy sufficient decrease" warnings without hitting the
    # solver iteration cap); tolerance tuning does not clear them, so they are silenced
    # with `verbosity = 0`.
    sol = integrate(idae, VSPARK(SPARKLobABC(3)); verbosity=0, warn_iterations=0)
    @test relative_maximum_error(sol.q, ref.q) < 4E-10

    sol = integrate(idae, VSPARK(SPARKLobABC(4)))
    @test relative_maximum_error(sol.q, ref.q) < 8E-15


    # order-reduced vs s=4; the solver stalls at a few steps regardless of tolerance
    # or solver choice, so its warnings are silenced with `verbosity = 0` rather than
    # tuned away.
    sol = integrate(idae, VSPARK(SPARKLobABD(3)); verbosity=0, warn_iterations=0)
    @test relative_maximum_error(sol.q, ref.q) < 8E-6

    # converges to 5.2E-12. Under SimpleSolvers 0.9.2 this stalled one step just above
    # machine precision and warned regardless of f_abstol. The SPARK method default is
    # now f_abstol = 8e-15 (see src/spark/abstract.jl), one order above this residual
    # floor, so under 0.10 the solve converges silently — no warning left to suppress.
    sol = integrate(idae, VSPARK(SPARKLobABD(4)))
    @test relative_maximum_error(sol.q, ref.q) < 8E-12


    sol = integrate(idae, VSPARK(SPARKLobattoIIIAIIIB(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, VSPARK(SPARKLobattoIIIAIIIB(4)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(idae, VSPARK(SPARKLobattoIIIBIIIA(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, VSPARK(SPARKLobattoIIIBIIIA(4)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    # SPARKLobattoIIIBIIIA(2) was formerly @test_broken (a singular stage system). Under
    # SimpleSolvers 0.10 its solve now converges, so this case passes and is asserted.
    sol = integrate(idae, VSPARK(SPARKLobattoIIIBIIIA(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # --- known-broken (see docs/src/audit.md) --- assert the mechanism per case.

    # Singular stage system at s=2 (the coinciding tableau pair) → SingularException.
    @test_throws SingularException integrate(idae, VSPARK(SPARKGLRK(1)))
    @test_throws SingularException integrate(idae, VSPARK(SPARKGLRK(2)))
    @test_throws SingularException integrate(idae, VSPARK(SPARKLobattoIIIAIIIB(2)))
    # SPARKLobABC(2) stalls at the residual floor, SPARKLobABD(2) converges to an
    # order-deficient answer (meas 0.06) — both marginal, so kept as @test_broken.
    @test_broken relative_maximum_error(integrate(idae, VSPARK(SPARKLobABC(2)); verbosity=0, warn_iterations=0).q, ref.q) < 1E-6
    @test_broken relative_maximum_error(integrate(idae, VSPARK(SPARKLobABD(2)); verbosity=0, warn_iterations=0).q, ref.q) < 1E-6


    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIAIIIB(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIAIIIB(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIAIIIB(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIBIIIA(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIBIIIA(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIBIIIA(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


@testset "$(rpad("VSPARK integrators with projection on primary constraint",80))" begin

    ## VSPARKprimary Integrators ###

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauVSPARKGLRKpSymplectic(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKGLRKpSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11


    sol = integrate(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(3))
    @test relative_maximum_error(sol.q, ref.q) < 8E-11

    sol = integrate(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

    sol = integrate(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    # --- known-broken (see docs/src/audit.md) ---

    # order reduction (meas 7.95E-7, no better than s=2); the solve stalls at the floor.
    @test nlsolve_outcome(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(3)).stalled

end


@testset "$(rpad("VSPARK integrators with projection on secondary constraint",80))" begin

    ## VSPARKsecondary Integrators ###

    sol = integrate(ldae, TableauVSPARKLobattoIIIAB(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIAB(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIAB(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIIBA(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIBA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIBA(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIICC̄(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIICC̄(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIICC̄(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIIC̄C(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIC̄C(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIC̄C(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIID(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIID(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIID(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-7

    sol = integrate(ldae, TableauVSPARKLobattoIIIE(3))
    @test relative_maximum_error(sol.q, ref.q) < 8E-12

    sol = integrate(ldae, TableauVSPARKLobattoIIIE(4))
    @test relative_maximum_error(sol.q, ref.q) < 8E-16


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIAB(1))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIAB(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-12

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIAB(3))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIBA(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIBA(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-12

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIBA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIICC̄(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIICC̄(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIICC̄(3))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIC̄C(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIC̄C(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIC̄C(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIID(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIID(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIID(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIE(1))
    @test relative_maximum_error(sol.q, ref.q) < 8E-7

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-12

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIE(3))
    @test relative_maximum_error(sol.q, ref.q) < 8E-16

end


@testset "$(rpad("HPARK integrators",80))" begin

    sol = integrate(pdae, TableauHPARKGLRK(1))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(pdae, TableauHPARKGLRK(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(pdae, TableauHPARKLobattoIIIAIIIB(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3
    # TODO: Check errors and large number of solver iterations !!!

    sol = integrate(pdae, TableauHPARKLobattoIIIBIIIA(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-2
    # TODO: Check errors and large number of solver iterations !!!


    # --- known-broken (see docs/src/audit.md) ---

    # The per-step solves converge cleanly, but the method is unstable on this problem:
    # the trajectory grows instead of tracking the reference (meas 20.0 / 11.3 / 1.61 / 0.15).
    @test nlsolve_outcome(pdae, TableauHPARKLobattoIIIAIIIB(2)).converged
    @test nlsolve_outcome(pdae, TableauHPARKLobattoIIIAIIIB(3)).converged
    @test nlsolve_outcome(pdae, TableauHPARKLobattoIIIBIIIA(2)).converged
    @test nlsolve_outcome(pdae, TableauHPARKLobattoIIIBIIIA(3)).converged

end


# TODO: Replace pdae with hdae !!!

@testset "$(rpad("HSPARK integrators",80))" begin

    # println("HSPARK(SPARKGLRK(1))")
    sol = integrate(pdae, HSPARK(SPARKGLRK(1)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    # println("HSPARK(SPARKGLRK(2))")
    sol = integrate(pdae, HSPARK(SPARKGLRK(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11


    # println("HSPARK(SPARKLobABC(2))")
    sol = integrate(pdae, HSPARK(SPARKLobABC(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # println("HSPARK(SPARKLobABC(3))")
    sol = integrate(pdae, HSPARK(SPARKLobABC(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    # println("HSPARK(SPARKLobABC(4))")
    sol = integrate(pdae, HSPARK(SPARKLobABC(4)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    # println("HSPARK(SPARKLobABD(2))")
    sol = integrate(pdae, HSPARK(SPARKLobABD(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # println("HSPARK(SPARKLobABD(3))")
    sol = integrate(pdae, HSPARK(SPARKLobABD(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    # println("HSPARK(SPARKLobABD(4))")
    sol = integrate(pdae, HSPARK(SPARKLobABD(4)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    # order reduction: accuracy plateaus ~2E-4, does not improve with s
    sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIAIIIB(1)))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIAIIIB(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIAIIIB(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-4


    sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIBIIIA(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIBIIIA(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIBIIIA(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-4


    # --- known-broken (see docs/src/audit.md) ---

    # Singular stage system (the coinciding Lobatto pair) → SingularException at every s.
    @test_throws SingularException integrate(pdae, HSPARK(SPARKLobattoIIIAIIIB(2)))
    @test_throws SingularException integrate(pdae, HSPARK(SPARKLobattoIIIAIIIB(3)))
    @test_throws SingularException integrate(pdae, HSPARK(SPARKLobattoIIIAIIIB(4)))
    @test_throws SingularException integrate(pdae, HSPARK(SPARKLobattoIIIBIIIA(2)))
    @test_throws SingularException integrate(pdae, HSPARK(SPARKLobattoIIIBIIIA(3)))
    @test_throws SingularException integrate(pdae, HSPARK(SPARKLobattoIIIBIIIA(4)))

end


@testset "$(rpad("HSPARK integrators with projection on primary constraint",80))" begin

    ### HSPARKprimary Integrators ###

    sol = integrate(pdae, TableauHSPARKGLRKpSymmetric(1))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKGLRKpSymmetric(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(3))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(3))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

end


@testset "$(rpad("HSPARK integrators with projection on secondary constraint",80))" begin

    ### HSPARKsecondary Integrators ###

    # HSPARKsecondary is EXPERIMENTAL and remains @test_broken (see the "SPARK
    # submodule" pass in docs/src/audit.md). Two implementation bugs were fixed:
    #   * the momentum projection coefficients a_p_2 / a_p_3 in getTableauHSPARK were
    #     s×s instead of s×σ, so the GLRK variants raised a BoundsError — now built as
    #     the conjugate-symplectic s×σ partners of α_q_2 / α_q_3;
    #   * the null-vector residual/component code was commented out while the cache
    #     still allocated the μ unknown, leaving an unconstrained (zero) Jacobian row —
    #     now re-enabled, matching the working VSPARKsecondary.
    # A residual singularity remains in the ω secondary-constraint block: every variant
    # still raises a SingularException, so all cases stay @test_broken. The solves now
    # iterate before failing (they no longer abort immediately), so they are silenced
    # with `verbosity = 0`.
    # The residual singularity in the ω secondary-constraint block raises a
    # SingularException for every variant and every s — that is the assertion.
    for s in (2, 3, 4)
        @test_throws SingularException integrate(hdae, TableauHSPARKLobattoIIIAB(s))
        @test_throws SingularException integrate(hdae, TableauHSPARKLobattoIIIBA(s))
        @test_throws SingularException integrate(hdae, TableauHSPARKLobattoIIID(s))
        @test_throws SingularException integrate(hdae, TableauHSPARKLobattoIIIE(s))

        @test_throws SingularException integrate(hdae, TableauHSPARKGLRKLobattoIIIAB(s))
        @test_throws SingularException integrate(hdae, TableauHSPARKGLRKLobattoIIIBA(s))
        @test_throws SingularException integrate(hdae, TableauHSPARKGLRKLobattoIIID(s))
        @test_throws SingularException integrate(hdae, TableauHSPARKGLRKLobattoIIIE(s))
    end

end

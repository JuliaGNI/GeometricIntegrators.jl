using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
import GeometricProblems.PointVortices as PointVortices
using Test

include("verification_utilities.jl")

# Formal Lagrangian Runge-Kutta methods on the (degenerate) Lotka-Volterra 2d
# Lagrangian system, referenced against a high-order solution of the equivalent ODE.
#
# There is NO order theorem for FLRK anywhere in the reference: the numerical-experiment
# sections are empty stubs and the backward-error analysis covers only s = 1. The
# expected orders below are therefore the order of the underlying Gauss tableau (2s),
# passed explicitly and *measured* rather than taken on trust.
#
# `v̄` must depend on q only for FLRK to be meaningful, which holds for
# `LotkaVolterra2d.lodeproblem` (`v̄ = v(t,q)`) but NOT for
# `HarmonicOscillator.lodeproblem` (`v̄ = p₁`, so ∂v̄/∂q ≡ 0).

const T = 1.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

build(Δt) = lodeproblem(q₀; timespan=(0.0, T), timestep=Δt, parameters=params)
ref(prob) = integrate(odeproblem(q₀; timespan=timespan(prob), timestep=timestep(prob), parameters=params), Gauss(8))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, r) = relative_maximum_error(sol.q, r.q)

# The adjoint variable should track ϑ(q): the momentum-modified formal Lagrangian is
# built precisely so that p = ϑ(q) solves the adjoint equations.
function emp(sol, r)
    maximum(maximum(abs, sol.p[i] .- LotkaVolterra2d.ϑ(sol.t[i], sol.q[i]))
            for i in eachindex(sol.q))
end

@testset "$(rpad("Formal Lagrangian Runge-Kutta convergence",80))" begin

    @testset "position" begin
        test_convergence_order(build, FLRK(Gauss(1)), steps(4, 4); reference=ref, errormetric=emq, expected=2, label="FLRK(Gauss(1)) q")
        test_convergence_order(build, FLRK(Gauss(2)), steps(4, 4); reference=ref, errormetric=emq, expected=4, label="FLRK(Gauss(2)) q")
        test_convergence_order(build, FLRK(Gauss(3)), steps(4, 3); reference=ref, errormetric=emq, expected=6, label="FLRK(Gauss(3)) q")
    end

    @testset "adjoint momentum" begin
        test_convergence_order(build, FLRK(Gauss(1)), steps(4, 4); reference=ref, errormetric=emp, expected=2, label="FLRK(Gauss(1)) p")
        test_convergence_order(build, FLRK(Gauss(2)), steps(4, 4); reference=ref, errormetric=emp, expected=4, label="FLRK(Gauss(2)) p")
        test_convergence_order(build, FLRK(Gauss(3)), steps(4, 3); reference=ref, errormetric=emp, expected=6, label="FLRK(Gauss(3)) p")
    end

    # The point vortices are the problem this method exists for: `GeometricProblems`
    # provides `lodeproblem_formal_lagrangian` specifically for the formal-Lagrangian
    # formulation, and its v̄ depends on q only, as FLRK requires. Measured orders here are
    # exactly 2s in both components — 2.02/4.00/6.00 in q and 2.03/4.00/6.01 in p.
    @testset "point vortices (formal Lagrangian formulation)" begin
        Tpv = 0.5
        pvbuild(Δt) = PointVortices.lodeproblem_formal_lagrangian(; timespan=(0.0, Tpv), timestep=Δt)
        pvref(prob) = integrate(PointVortices.odeproblem(; timespan=timespan(prob), timestep=timestep(prob)), Gauss(8))
        pvsteps(n0, k) = Tpv ./ (n0 .* 2 .^ (0:k))
        pvemp(sol, r) = maximum(maximum(abs, sol.p[i] .- PointVortices.ϑ(sol.q[i])) for i in eachindex(sol.q))

        for (s, expected) in ((1, 2), (2, 4), (3, 6))
            test_convergence_order(pvbuild, FLRK(Gauss(s)), pvsteps(16, 3); reference=pvref, errormetric=emq, expected=expected, label="FLRK(Gauss($s)) q [vortices]")
            test_convergence_order(pvbuild, FLRK(Gauss(s)), pvsteps(16, 3); reference=pvref, errormetric=pvemp, expected=expected, minpoints=2, label="FLRK(Gauss($s)) p [vortices]")
        end
    end

    # The position phase is an implicit Runge-Kutta method applied to q̇ = v̄(q), so it
    # must agree with `Gauss(s)` on the ODE form to round-off at every timestep. This is
    # a far stronger statement than any error bound: it pins the whole q-phase.
    #
    # Not bit-for-bit: FLRK seeds the Newton iteration from `v̄(t,q,p)` on top of the
    # Hermite extrapolation, so the converged iterates can differ in the last bit. The
    # measured discrepancy is exactly 1 ulp (2.2E-16 absolute, 1.6E-16 relative) at the
    # coarsest step and identically zero at the finer ones.
    @testset "consistency with the underlying Runge-Kutta method" begin
        for Δt in steps(4, 2), s in 1:3
            prob = build(Δt)
            odeprob = odeproblem(q₀; timespan=timespan(prob), timestep=timestep(prob), parameters=params)
            @test relative_maximum_error(integrate(prob, FLRK(Gauss(s))).q,
                                        integrate(odeprob, Gauss(s)).q) < 4E-16
        end
    end

end

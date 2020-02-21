
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators: initial_guess!, jacobian!, update_params!
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using GeometricIntegrators.Utils
using Test

set_config(:nls_stol_break, 1E3)

using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem: Δt, nt, refx, refq, refp, k

ode  = harmonic_oscillator_ode()
pode = harmonic_oscillator_pode()
sode = harmonic_oscillator_sode()
iode = harmonic_oscillator_iode()
idae = harmonic_oscillator_idae()
pdae = harmonic_oscillator_pdae()


@testset "$(rpad("Runge-Kutta integrators",80))" begin

    @test typeof(Integrator(ode, getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
    @test typeof(Integrator(ode, getTableauCrouzeix(), Δt)) <: IntegratorDIRK
    @test typeof(Integrator(ode, getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK
    @test typeof(Integrator(ode, getTableauGLRK(1), Δt)) <: IntegratorFIRK

end


@testset "$(rpad("Explicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, getTableauExplicitEuler(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-2

    sol = integrate(ode, getTableauExplicitMidpoint(), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-3

    sol = integrate(ode, getTableauERK4(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-7

end


function test_firk_jacobian(ode, tableau, Δt, nt; atol=eps())
    int1 = IntegratorFIRK(ode, tableau, Δt; exact_jacobian=false)
    int2 = IntegratorFIRK(ode, tableau, Δt; exact_jacobian=true)

    sol = Solution(ode, Δt, nt)

    asol1 = AtomicSolution(ode)
    asol2 = AtomicSolution(ode)

    get_initial_conditions!(sol, asol1, 1, 1)
    get_initial_conditions!(sol, asol2, 1, 1)

    update_params!(int1, asol1)
    update_params!(int2, asol2)

    initial_guess!(int1, asol1)
    initial_guess!(int2, asol2)

    computeJacobian(int1.solver.x, int1.solver.J, int1.solver.Jparams)
    computeJacobian(int2.solver.x, int2.solver.J, int2.solver.Jparams)

    @test int1.solver.J ≈ int2.solver.J atol=atol
end


@testset "$(rpad("Implicit Runge-Kutta integrators",80))" begin

    # check Jacobian of FIRK integrators

    int = Integrator(ode, getTableauImplicitMidpoint(), Δt)
    jacobian!(int.solver.x, int.solver.J, IntegratorCache(int.params), int.params)
    @test int.solver.J == - [1.0  0.0; 0.0  1.0] + Δt / 2 * [0.0  1.0; -k  0.0]

    test_firk_jacobian(ode, getTableauGLRK(1), Δt, nt)
    test_firk_jacobian(ode, getTableauGLRK(2), Δt, nt)
    test_firk_jacobian(ode, getTableauGLRK(3), Δt, nt)
    test_firk_jacobian(ode, getTableauGLRK(4), Δt, nt)


    # check integrators

    sol = integrate(ode, getTableauCrouzeix(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-5

    sol = integrate(ode, getTableauImplicitEuler(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-2

    sol = integrate(ode, getTableauImplicitMidpoint(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-4

    sol = integrate(ode, getTableauGLRK(1), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-4

    sol = integrate(ode, getTableauGLRK(2), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-7

    sol = integrate(ode, getTableauGLRK(3), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-11

    sol = integrate(ode, getTableauGLRK(4), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-15

    sol = integrate(ode, getTableauGLRK(5), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-15

    sol = integrate(ode, getTableauGLRK(6), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-15

    sol = integrate(ode, getTableauGLRK(7), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-15

    sol = integrate(ode, getTableauSRK3(), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-7


    # Test integration with exact Jacobian

    int = IntegratorFIRK(ode, getTableauGLRK(1), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 5E-4

    int = IntegratorFIRK(ode, getTableauGLRK(2), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 1E-7

    int = IntegratorFIRK(ode, getTableauGLRK(3), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorFIRK(ode, getTableauGLRK(4), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 1E-15

end


@testset "$(rpad("Partitioned Runge-Kutta integrators",80))" begin

    psol = integrate(pode, getTableauSymplecticEulerA(), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-2
    @test rel_err(psol.p, refp) < 1E-3

    psol = integrate(pode, getTableauSymplecticEulerB(), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-2
    @test rel_err(psol.p, refp) < 1E-3

    psol = integrate(pode, TableauEPRK(:prk4, 4, getTableauERK4().q), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-7
    @test rel_err(psol.p, refp) < 5E-7

    psol = integrate(pode, TableauIPRK(:pglrk, 2, getCoefficientsGLRK(1)), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-4
    @test rel_err(psol.p, refp) < 5E-4

    psol = integrate(pode, TableauIPRK(:pglrk, 4, getCoefficientsGLRK(2)), Δt, nt)
    @test rel_err(psol.q, refq) < 1E-7
    @test rel_err(psol.p, refp) < 1E-7

    psol = integrate(pode, TableauIPRK(:pglrk, 6, getCoefficientsGLRK(3)), Δt, nt)
    @test rel_err(psol.q, refq) < 1E-11
    @test rel_err(psol.p, refp) < 1E-11

end


@testset "$(rpad("Special Runge-Kutta integrators",80))" begin

    pgint = IntegratorPGLRK(ode, getCoefficientsPGLRK(2), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 4E-8

    pgint = IntegratorPGLRK(ode, getCoefficientsPGLRK(3), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 2E-12

    pgint = IntegratorPGLRK(ode, getCoefficientsPGLRK(4), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 8E-16

end

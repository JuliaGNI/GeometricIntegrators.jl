
using GeometricBase.Config
using GeometricBase.Utils
using GeometricEquations
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators: initial_guess!, jacobian!, update_params!
using GeometricIntegrators.Solutions
using GeometricIntegrators.Tableaus
using GeometricProblems.HarmonicOscillator
using SimpleSolvers
using Test

using GeometricProblems.HarmonicOscillator: Δt, nt, refx, refq, refp, k

ode  = harmonic_oscillator_ode()
pode = harmonic_oscillator_pode()
sode = harmonic_oscillator_sode()
iode = harmonic_oscillator_iode()
idae = harmonic_oscillator_idae()
hode = harmonic_oscillator_hode()
pdae = harmonic_oscillator_pdae()



@testset "$(rpad("Runge-Kutta integrators",80))" begin

    @test typeof(Integrator(ode, TableauExplicitMidpoint(), Δt)) <: IntegratorERK
    @test typeof(Integrator(ode, TableauCrouzeix(), Δt)) <: IntegratorDIRK
    @test typeof(Integrator(ode, TableauImplicitMidpoint(), Δt)) <: IntegratorFIRK
    @test typeof(Integrator(pode, PartitionedTableau(:eprk_midpoint, TableauExplicitMidpoint()), Δt)) <: IntegratorEPRK
    @test typeof(Integrator(pode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()), Δt)) <: IntegratorIPRK
    @test typeof(Integrator(iode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()), Δt)) <: IntegratorPRKimplicit

    @test_logs (:warn, r"Initializing IntegratorDIRK with explicit tableau heun.*") IntegratorDIRK(ode, TableauHeun2(), Δt)
    @test_logs (:warn, r"Initializing IntegratorFIRK with explicit tableau heun.*") IntegratorFIRK(ode, TableauHeun2(), Δt)
    @test_logs (:warn, r"Initializing IntegratorFIRK with diagonally implicit tableau crouzeix.*") IntegratorFIRK(ode, TableauCrouzeix(), Δt)

    int_erk  = Integrator(ode, TableauExplicitMidpoint(), Δt)
    int_dirk = Integrator(ode, TableauCrouzeix(), Δt)
    int_firk = Integrator(ode, TableauImplicitMidpoint(), Δt)
    int_eprk = Integrator(pode, PartitionedTableau(:eprk_midpoint, TableauExplicitMidpoint()), Δt)
    int_iprk = Integrator(pode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()), Δt)

    @test ndims(ode) == ndims(int_erk)
    @test ndims(ode) == ndims(int_dirk)
    @test ndims(ode) == ndims(int_firk)
    @test ndims(pode) == ndims(int_eprk)
    @test ndims(pode) == ndims(int_iprk)

end


@testset "$(rpad("Explicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, TableauExplicitEuler(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-2

    sol = integrate(ode, TableauExplicitMidpoint(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-3

    sol = integrate(ode, TableauRK4(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-7

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

    int = Integrator(ode, TableauImplicitMidpoint(), Δt)
    jacobian!(int.solver.x, int.solver.J, IntegratorCache(int.params), int.params)
    @test int.solver.J == - [1.0  0.0; 0.0  1.0] + Δt / 2 * [0.0  1.0; -k  0.0]

    test_firk_jacobian(ode, TableauGauss(1), Δt, nt)
    test_firk_jacobian(ode, TableauGauss(2), Δt, nt)
    test_firk_jacobian(ode, TableauGauss(3), Δt, nt)
    test_firk_jacobian(ode, TableauGauss(4), Δt, nt)


    # check integrators

    sol = integrate(ode, TableauCrouzeix(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-5

    sol = integrate(ode, TableauImplicitEuler(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-2

    sol = integrate(ode, TableauImplicitMidpoint(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-4

    sol = integrate(ode, TableauGauss(1), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-4

    sol = integrate(ode, TableauGauss(2), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-7

    sol = integrate(ode, TableauGauss(3), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    sol = integrate(ode, TableauGauss(4), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(5), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(6), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(7), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(8), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-15

    sol = integrate(ode, TableauSRK3(), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-7


    # Test integration with exact Jacobian

    int = IntegratorFIRK(ode, TableauGauss(1), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-4

    int = IntegratorFIRK(ode, TableauGauss(2), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-7

    int = IntegratorFIRK(ode, TableauGauss(3), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorFIRK(ode, TableauGauss(4), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

end


@testset "$(rpad("Partitioned Runge-Kutta integrators",80))" begin

    # psol = integrate(pode, TableauSymplecticEulerA(), Δt, nt)
    psol = integrate(pode, PartitionedTableau(:seulerA, TableauExplicitEuler(), TableauImplicitEuler()), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 5E-2
    @test relative_maximum_error(psol.p, refp) < 1E-3

    # psol = integrate(pode, TableauSymplecticEulerB(), Δt, nt)
    psol = integrate(pode, PartitionedTableau(:seulerB, TableauImplicitEuler(), TableauExplicitEuler()), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 5E-2
    @test relative_maximum_error(psol.p, refp) < 1E-3

    psol = integrate(pode, PartitionedTableau(:prk4, TableauRK4()), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 5E-7
    @test relative_maximum_error(psol.p, refp) < 5E-7

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(1)), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 5E-4
    @test relative_maximum_error(psol.p, refp) < 5E-4

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(2)), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 1E-7
    @test relative_maximum_error(psol.p, refp) < 1E-7

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(3)), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 1E-11
    @test relative_maximum_error(psol.p, refp) < 1E-11

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(4)), Δt, nt)
    @test relative_maximum_error(psol.q, refq) < 1E-15
    @test relative_maximum_error(psol.p, refp) < 1E-15

end


@testset "$(rpad("Special Runge-Kutta integrators",80))" begin

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(2), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test relative_maximum_error(pgsol.q, refx) < 6E-6

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(3), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test relative_maximum_error(pgsol.q, refx) < 2E-12

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(4), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test relative_maximum_error(pgsol.q, refx) < 8E-16

end


@testset "$(rpad("Integrate PODE and HODE with ODE Runge-Kutta integrators",80))" begin
    for s in 1:4
        code = convert(ODE, pode)
        csol = integrate(code, TableauGauss(s), Δt, nt)
        psol = integrate(pode, PartitionedTableauGauss(s), Δt, nt)

        @test csol.q[1, end] == psol.q[1, end]
        @test csol.q[2, end] == psol.p[1, end]
    end

    for s in 1:4
        code = convert(ODE, hode)
        csol = integrate(code, TableauGauss(s), Δt, nt)
        hsol = integrate(hode, PartitionedTableauGauss(s), Δt, nt)

        @test csol.q[1, end] == hsol.q[1, end]
        @test csol.q[2, end] == hsol.p[1, end]
    end

    for s in 1:4
        code = convert(PODE, hode)
        csol = integrate(code, PartitionedTableauGauss(s), Δt, nt)
        hsol = integrate(hode, PartitionedTableauGauss(s), Δt, nt)

        @test csol.q[1, end] == hsol.q[1, end]
        @test csol.p[1, end] == hsol.p[1, end]
    end
end

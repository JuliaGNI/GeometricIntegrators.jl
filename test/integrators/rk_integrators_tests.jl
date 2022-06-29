
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

    @test typeof(Integrator(ode, TableauExplicitMidpoint())) <: IntegratorERK
    @test typeof(Integrator(ode, TableauCrouzeix())) <: IntegratorDIRK
    @test typeof(Integrator(ode, TableauImplicitMidpoint())) <: IntegratorFIRK
    @test typeof(Integrator(pode, PartitionedTableau(:eprk_midpoint, TableauExplicitMidpoint()))) <: IntegratorEPRK
    @test typeof(Integrator(pode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()))) <: IntegratorIPRK
    @test typeof(Integrator(iode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()))) <: IntegratorPRKimplicit

    @test_logs (:warn, r"Initializing IntegratorDIRK with explicit tableau heun.*") IntegratorDIRK(ode, TableauHeun2())
    @test_logs (:warn, r"Initializing IntegratorFIRK with explicit tableau heun.*") IntegratorFIRK(ode, TableauHeun2())
    @test_logs (:warn, r"Initializing IntegratorFIRK with diagonally implicit tableau crouzeix.*") IntegratorFIRK(ode, TableauCrouzeix())

    int_erk  = Integrator(ode, TableauExplicitMidpoint())
    int_dirk = Integrator(ode, TableauCrouzeix())
    int_firk = Integrator(ode, TableauImplicitMidpoint())
    int_eprk = Integrator(pode, PartitionedTableau(:eprk_midpoint, TableauExplicitMidpoint()))
    int_iprk = Integrator(pode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()))

    @test ndims(ode) == ndims(int_erk)
    @test ndims(ode) == ndims(int_dirk)
    @test ndims(ode) == ndims(int_firk)
    @test ndims(pode) == ndims(int_eprk)
    @test ndims(pode) == ndims(int_iprk)

end


@testset "$(rpad("Explicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, TableauExplicitEuler(), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-2

    sol = integrate(ode, TableauExplicitMidpoint(), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-3

    sol = integrate(ode, TableauRK4(), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-7

end


function test_firk_jacobian(ode, tableau, nt; atol=eps())
    int1 = IntegratorFIRK(ode, tableau; exact_jacobian=false)
    int2 = IntegratorFIRK(ode, tableau; exact_jacobian=true)

    sol = Solution(ode, nt)

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

    int = Integrator(ode, TableauImplicitMidpoint())
    jacobian!(int.solver.x, int.solver.J, IntegratorCache(int.params), int.params)
    @test int.solver.J == - [1.0  0.0; 0.0  1.0] + Δt / 2 * [0.0  1.0; -k  0.0]

    test_firk_jacobian(ode, TableauGauss(1), nt)
    test_firk_jacobian(ode, TableauGauss(2), nt)
    test_firk_jacobian(ode, TableauGauss(3), nt)
    test_firk_jacobian(ode, TableauGauss(4), nt)


    # check integrators

    sol = integrate(ode, TableauCrouzeix(), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-5

    sol = integrate(ode, TableauImplicitEuler(), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-2

    sol = integrate(ode, TableauImplicitMidpoint(), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-4

    sol = integrate(ode, TableauGauss(1), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-4

    sol = integrate(ode, TableauGauss(2), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-7

    sol = integrate(ode, TableauGauss(3), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    sol = integrate(ode, TableauGauss(4), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(5), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(6), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(7), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGauss(8), nt)
    @test relative_maximum_error(sol.q, refx) < 5E-15

    sol = integrate(ode, TableauSRK3(), nt)
    @test relative_maximum_error(sol.q, refx) < 1E-7


    # Test integration with exact Jacobian

    int = IntegratorFIRK(ode, TableauGauss(1); exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 5E-4

    int = IntegratorFIRK(ode, TableauGauss(2); exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-7

    int = IntegratorFIRK(ode, TableauGauss(3); exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorFIRK(ode, TableauGauss(4); exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

end


@testset "$(rpad("Partitioned Runge-Kutta integrators",80))" begin

    # psol = integrate(pode, TableauSymplecticEulerA(), nt)
    psol = integrate(pode, PartitionedTableau(:seulerA, TableauExplicitEuler(), TableauImplicitEuler()), nt)
    @test relative_maximum_error(psol.q, refq) < 5E-2
    @test relative_maximum_error(psol.p, refp) < 1E-3

    # psol = integrate(pode, TableauSymplecticEulerB(), nt)
    psol = integrate(pode, PartitionedTableau(:seulerB, TableauImplicitEuler(), TableauExplicitEuler()), nt)
    @test relative_maximum_error(psol.q, refq) < 5E-2
    @test relative_maximum_error(psol.p, refp) < 1E-3

    psol = integrate(pode, PartitionedTableau(:prk4, TableauRK4()), nt)
    @test relative_maximum_error(psol.q, refq) < 5E-7
    @test relative_maximum_error(psol.p, refp) < 5E-7

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(1)), nt)
    @test relative_maximum_error(psol.q, refq) < 5E-4
    @test relative_maximum_error(psol.p, refp) < 5E-4

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(2)), nt)
    @test relative_maximum_error(psol.q, refq) < 1E-7
    @test relative_maximum_error(psol.p, refp) < 1E-7

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(3)), nt)
    @test relative_maximum_error(psol.q, refq) < 1E-11
    @test relative_maximum_error(psol.p, refp) < 1E-11

    psol = integrate(pode, PartitionedTableau(:pglrk, TableauGauss(4)), nt)
    @test relative_maximum_error(psol.q, refq) < 1E-15
    @test relative_maximum_error(psol.p, refp) < 1E-15

end


@testset "$(rpad("Special Runge-Kutta integrators",80))" begin

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(2))
    pgsol = integrate(ode, pgint, nt)
    @test relative_maximum_error(pgsol.q, refx) < 6E-6

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(3))
    pgsol = integrate(ode, pgint, nt)
    @test relative_maximum_error(pgsol.q, refx) < 2E-12

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(4))
    pgsol = integrate(ode, pgint, nt)
    @test relative_maximum_error(pgsol.q, refx) < 8E-16

end


@testset "$(rpad("Integrate PODE and HODE with ODE Runge-Kutta integrators",80))" begin
    for s in 1:4
        code = convert(ODEProblem, pode)
        csol = integrate(code, TableauGauss(s), nt)
        psol = integrate(pode, PartitionedTableauGauss(s), nt)

        @test csol.q[1, end] == psol.q[1, end]
        @test csol.q[2, end] == psol.p[1, end]
    end

    for s in 1:4
        code = convert(ODEProblem, hode)
        csol = integrate(code, TableauGauss(s), nt)
        hsol = integrate(hode, PartitionedTableauGauss(s), nt)

        @test csol.q[1, end] == hsol.q[1, end]
        @test csol.q[2, end] == hsol.p[1, end]
    end

    for s in 1:4
        code = convert(PODEProblem, hode)
        csol = integrate(code, PartitionedTableauGauss(s), nt)
        hsol = integrate(hode, PartitionedTableauGauss(s), nt)

        @test csol.q[1, end] == hsol.q[1, end]
        @test csol.p[1, end] == hsol.p[1, end]
    end
end

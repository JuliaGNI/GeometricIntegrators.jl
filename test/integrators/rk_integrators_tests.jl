using GeometricIntegrators
# using GeometricIntegrators.Integrators: initial_guess!, jacobian!, update_params!
using Test

using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: reference_solution, reference_solution_q, reference_solution_p

ode  = odeproblem()
pode = podeproblem()
sode = sodeproblem()
iode = iodeproblem()
idae = idaeproblem()
hode = hodeproblem()
pdae = pdaeproblem()

ref  = exact_solution(ode)
pref = exact_solution(pode)

# @testset "$(rpad("Runge-Kutta integrators",80))" begin

#     @test typeof(AbstractIntegrator(ode, TableauExplicitMidpoint())) <: IntegratorERK
#     @test typeof(AbstractIntegrator(ode, TableauCrouzeix())) <: IntegratorDIRK
#     @test typeof(AbstractIntegrator(ode, TableauImplicitMidpoint())) <: IntegratorIRK
#     @test typeof(AbstractIntegrator(pode, PartitionedTableau(:eprk_midpoint, TableauExplicitMidpoint()))) <: IntegratorEPRK
#     @test typeof(AbstractIntegrator(pode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()))) <: IntegratorIPRK
#     @test typeof(AbstractIntegrator(iode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()))) <: IntegratorPRKimplicit

#     @test_logs (:warn, r"Initializing IntegratorDIRK with explicit tableau Heun2.*") IntegratorDIRK(ode, TableauHeun2())
#     @test_logs (:warn, r"Initializing IntegratorIRK with explicit tableau Heun2.*") IntegratorIRK(ode, TableauHeun2())
#     @test_logs (:warn, r"Initializing IntegratorIRK with diagonally implicit tableau Crouzeix.*") IntegratorIRK(ode, TableauCrouzeix())

#     int_erk  = AbstractIntegrator(ode, TableauExplicitMidpoint())
#     int_dirk = AbstractIntegrator(ode, TableauCrouzeix())
#     int_irk = AbstractIntegrator(ode, TableauImplicitMidpoint())
#     int_eprk = AbstractIntegrator(pode, PartitionedTableau(:eprk_midpoint, TableauExplicitMidpoint()))
#     int_iprk = AbstractIntegrator(pode, PartitionedTableau(:iprk_midpoint, TableauImplicitMidpoint()))

#     @test ndims(ode) == ndims(int_erk)
#     @test ndims(ode) == ndims(int_dirk)
#     @test ndims(ode) == ndims(int_irk)
#     @test ndims(pode) == ndims(int_eprk)
#     @test ndims(pode) == ndims(int_iprk)

# end


@testset "$(rpad("Explicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, ExplicitEulerRK())
    @test relative_maximum_error(sol, ref).q < 5E-2

    sol = integrate(ode, ExplicitMidpoint())
    @test relative_maximum_error(sol, ref).q < 1E-3

    sol = integrate(ode, RK4())
    @test relative_maximum_error(sol, ref).q < 5E-7

end


@testset "$(rpad("Implicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, ImplicitEulerRK())
    @test relative_maximum_error(sol, ref).q < 5E-2

    sol = integrate(ode, ImplicitMidpoint())
    @test relative_maximum_error(sol, ref).q < 5E-4

    sol = integrate(ode, SRK3())
    @test relative_maximum_error(sol, ref).q < 1E-7

    sol = integrate(ode, Gauss(1))
    @test relative_maximum_error(sol, ref).q < 5E-4

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol, ref).q < 1E-7

    sol = integrate(ode, Gauss(3))
    @test relative_maximum_error(sol, ref).q < 1E-11

    sol = integrate(ode, Gauss(4))
    @test relative_maximum_error(sol, ref).q < 1E-15

    sol = integrate(ode, Gauss(5))
    @test relative_maximum_error(sol, ref).q < 1E-15

    sol = integrate(ode, Gauss(6))
    @test relative_maximum_error(sol, ref).q < 1E-15

    sol = integrate(ode, Gauss(7))
    @test relative_maximum_error(sol, ref).q < 1E-15

    sol = integrate(ode, Gauss(8))
    @test relative_maximum_error(sol, ref).q < 5E-15

end


# @testset "$(rpad("Implicit Runge-Kutta integrators with Block Jacobian",80))" begin

#     # Check Block Jacobian for IRK integrators

#     function test_irk_jacobian(ode, tableau; atol=eps())
#         int1 = IntegratorIRK(ode, tableau; exact_jacobian=false)
#         int2 = IntegratorIRK(ode, tableau; exact_jacobian=true)

#         solstep1 = SolutionStep(ode)
#         solstep2 = SolutionStep(ode)

#         initialize!(int1, solstep2)
#         initialize!(int1, solstep2)

#         update_params!(int1, solstep1)
#         update_params!(int2, solstep2)

#         initial_guess!(int1, solstep1)
#         initial_guess!(int2, solstep2)

#         computeJacobian(int1.solver.x, int1.solver.J, int1.solver.Jparams)
#         computeJacobian(int2.solver.x, int2.solver.J, int2.solver.Jparams)

#         @test int1.solver.J ≈ int2.solver.J atol=atol
#     end

#     int = AbstractIntegrator(ode, TableauImplicitMidpoint())
#     jacobian!(int.solver.x, int.solver.J, IntegratorCache(int.params), int.params)
#     @test int.solver.J == - [1.0  0.0; 0.0  1.0] + timestep(ode) / 2 * [0.0  1.0; -parameters(ode).k  0.0]

#     test_IRK_jacobian(ode, TableauGauss(1))
#     test_IRK_jacobian(ode, TableauGauss(2))
#     test_IRK_jacobian(ode, TableauGauss(3))
#     test_IRK_jacobian(ode, TableauGauss(4))


#     # Test integration with block Jacobian

#     int = IntegratorIRK(ode, TableauGauss(1); exact_jacobian=true)
#     sol = integrate(ode, int)
#     @test relative_maximum_error(sol, ref).q < 5E-4

#     int = IntegratorIRK(ode, TableauGauss(2); exact_jacobian=true)
#     sol = integrate(ode, int)
#     @test relative_maximum_error(sol, ref).q < 1E-7

#     int = IntegratorIRK(ode, TableauGauss(3); exact_jacobian=true)
#     sol = integrate(ode, int)
#     @test relative_maximum_error(sol, ref).q < 1E-11

#     int = IntegratorIRK(ode, TableauGauss(4); exact_jacobian=true)
#     sol = integrate(ode, int)
#     @test relative_maximum_error(sol, ref).q < 1E-15

# end


@testset "$(rpad("Diagonally Implicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, Crouzeix())
    @test relative_maximum_error(sol, ref).q < 5E-5

    sol = integrate(ode, CrankNicolson())
    @test relative_maximum_error(sol, ref).q < 5E-4

    sol = integrate(ode, KraaijevangerSpijker())
    @test relative_maximum_error(sol, ref).q < 1E-1

    sol = integrate(ode, QinZhang())
    @test relative_maximum_error(sol, ref).q < 1E-4

end


@testset "$(rpad("Explicit Partitioned Runge-Kutta integrators",80))" begin

    psol = integrate(pode, SymplecticEulerA())
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 5E-2
    # @test perr.p < 1E-3

    psol = integrate(pode, SymplecticEulerB())
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 5E-2
    # @test perr.p < 1E-3

    psol = integrate(pode, LobattoIIIAIIIB(2))
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 2E-4
    # @test perr.p < 5E-4

    psol = integrate(pode, LobattoIIIBIIIA(2))
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 2E-4
    # @test perr.p < 1E-3

    psol = integrate(pode, RK4())
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 2E-7
    # @test perr.p < 2E-7

end

@testset "$(rpad("Implicit Partitioned Runge-Kutta integrators",80))" begin

    psol = integrate(pode, Gauss(1))
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 5E-4
    # @test perr.p < 5E-4
    @test psol.q == integrate(pode, PartitionedGauss(1)).q
    @test psol.q == integrate(pode, ImplicitMidpoint()).q

    psol = integrate(pode, Gauss(2))
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 1E-7
    # @test perr.p < 1E-7
    @test psol.q == integrate(pode, PartitionedGauss(2)).q

    psol = integrate(pode, Gauss(3))
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 1E-11
    # @test perr.p < 1E-11
    @test psol.q == integrate(pode, PartitionedGauss(3)).q

    psol = integrate(pode, Gauss(4))
    perr = relative_maximum_error(psol, pref)
    @test perr.q < 1E-15
    # @test perr.p < 1E-15
    @test psol.q == integrate(pode, PartitionedGauss(4)).q

end


# @testset "$(rpad("Projected Gauss-Legendre Runge-Kutta integrators",80))" begin

#     pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(2))
#     pgsol = integrate(ode, pgint)
#     @test relative_maximum_error(pgsol.q, reference_solution) < 6E-6

#     pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(3))
#     pgsol = integrate(ode, pgint)
#     @test relative_maximum_error(pgsol.q, reference_solution) < 2E-12

#     pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(4))
#     pgsol = integrate(ode, pgint)
#     @test relative_maximum_error(pgsol.q, reference_solution) < 8E-16

# end


@testset "$(rpad("Integrate PODE and HODE with ODE Runge-Kutta integrators",80))" begin
    for s in 1:4
        code = convert(ODEProblem, pode)
        csol = integrate(code, Gauss(s))
        psol = integrate(pode, Gauss(s))

        @test csol.q[end][1] == psol.q[end][1]
        @test csol.q[end][2] == psol.p[end][1]
    end

    for s in 1:4
        code = convert(ODEProblem, hode)
        csol = integrate(code, Gauss(s))
        hsol = integrate(hode, Gauss(s))

        @test csol.q[end][1] == hsol.q[end][1]
        @test csol.q[end][2] == hsol.p[end][1]
    end

    for s in 1:4
        code = convert(PODEProblem, hode)
        csol = integrate(code, Gauss(s))
        hsol = integrate(hode, Gauss(s))

        @test csol.q[end][1] == hsol.q[end][1]
        @test csol.p[end][1] == hsol.p[end][1]
    end
end

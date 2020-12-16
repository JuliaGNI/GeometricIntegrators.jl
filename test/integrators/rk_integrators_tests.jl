
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators: initial_guess!, jacobian!, update_params!
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator
using Test

set_config(:nls_stol_break, 1E3)

using GeometricProblems.HarmonicOscillator: Δt, nt, refx, refq, refp, k

ode  = harmonic_oscillator_ode()
pode = harmonic_oscillator_pode()
sode = harmonic_oscillator_sode()
iode = harmonic_oscillator_iode()
idae = harmonic_oscillator_idae()
pdae = harmonic_oscillator_pdae()



@testset "$(rpad("Runge-Kutta coefficients",80))" begin

    s = 2

    g = CoefficientsGLRK(s)
    g̃ = get_symplectic_conjugate_coefficients(g)

    @test g.a ≈ g̃.a atol=1E-15
    @test g.b == g̃.b
    @test g.c == g̃.c

    @test compute_symplecticity_error(g) == zeros(s,s)
    @test check_symplecticity(g) == Array{Bool}(ones(s,s))
    @test check_symmetry(g) == Array{Bool}(ones(s,s))

    @test check_order_conditions_B(g,1) == true
    @test check_order_conditions_B(g,2) == true
    @test check_order_conditions_C(g,1) == Array{Bool}(ones(s))
    @test check_order_conditions_C(g,2) == Array{Bool}(ones(s))
    @test check_order_conditions_D(g,1) == Array{Bool}(ones(s))
    @test check_order_conditions_D(g,2) == Array{Bool}(ones(s))

    A = CoefficientsLobattoIIIA(s)
    B = CoefficientsLobattoIIIB(s)
    E = CoefficientsLobattoIIIE(s)
    Ã = get_symplectic_conjugate_coefficients(A)
    B̃ = get_symplectic_conjugate_coefficients(B)

    @test A.a ≈ B̃.a atol=1E-15
    @test A.b == Ã.b
    @test A.c == Ã.c

    @test B.a ≈ Ã.a atol=1E-15
    @test B.b == Ã.b
    @test B.c == Ã.c

    Â = symplecticize(A)
    B̂ = symplecticize(B)

    @test Â.a ≈ E.a atol=1E-15
    @test Â.b == E.b
    @test Â.c == E.c

    @test B̂.a ≈ E.a atol=1E-15
    @test B̂.b == E.b
    @test B̂.c == E.c

end


@testset "$(rpad("Runge-Kutta integrators",80))" begin

    @test typeof(Integrator(ode, TableauExplicitMidpoint(), Δt)) <: IntegratorERK
    @test typeof(Integrator(ode, TableauCrouzeix(), Δt)) <: IntegratorDIRK
    @test typeof(Integrator(ode, TableauImplicitMidpoint(), Δt)) <: IntegratorFIRK
    @test typeof(Integrator(pode, TableauEPRK(:eprk_midpoint, 2, TableauExplicitMidpoint().q), Δt)) <: IntegratorEPRK
    @test typeof(Integrator(pode, TableauIPRK(:iprk_midpoint, 2, TableauImplicitMidpoint().q), Δt)) <: IntegratorIPRK

    int_erk  = Integrator(ode, TableauExplicitMidpoint(), Δt)
    int_dirk = Integrator(ode, TableauCrouzeix(), Δt)
    int_firk = Integrator(ode, TableauImplicitMidpoint(), Δt)
    int_eprk = Integrator(pode, TableauEPRK(:eprk_midpoint, 2, TableauExplicitMidpoint().q), Δt)
    int_iprk = Integrator(pode, TableauIPRK(:iprk_midpoint, 2, TableauImplicitMidpoint().q), Δt)

    @test ndims(ode) == ndims(int_erk)
    @test ndims(ode) == ndims(int_dirk)
    @test ndims(ode) == ndims(int_firk)
    @test ndims(pode) == ndims(int_eprk)
    @test ndims(pode) == ndims(int_iprk)

end


@testset "$(rpad("Explicit Runge-Kutta integrators",80))" begin

    sol = integrate(ode, TableauExplicitEuler(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-2

    sol = integrate(ode, TableauExplicitMidpoint(), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-3

    sol = integrate(ode, TableauRK4(), Δt, nt)
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

    int = Integrator(ode, TableauImplicitMidpoint(), Δt)
    jacobian!(int.solver.x, int.solver.J, IntegratorCache(int.params), int.params)
    @test int.solver.J == - [1.0  0.0; 0.0  1.0] + Δt / 2 * [0.0  1.0; -k  0.0]

    test_firk_jacobian(ode, TableauGLRK(1), Δt, nt)
    test_firk_jacobian(ode, TableauGLRK(2), Δt, nt)
    test_firk_jacobian(ode, TableauGLRK(3), Δt, nt)
    test_firk_jacobian(ode, TableauGLRK(4), Δt, nt)


    # check integrators

    sol = integrate(ode, TableauCrouzeix(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-5

    sol = integrate(ode, TableauImplicitEuler(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-2

    sol = integrate(ode, TableauImplicitMidpoint(), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-4

    sol = integrate(ode, TableauGLRK(1), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-4

    sol = integrate(ode, TableauGLRK(2), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-7

    sol = integrate(ode, TableauGLRK(3), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-11

    sol = integrate(ode, TableauGLRK(4), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGLRK(5), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGLRK(6), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-15

    sol = integrate(ode, TableauGLRK(7), Δt, nt)
    @test rel_err(sol.q, refx) < 5E-15

    sol = integrate(ode, TableauSRK3(), Δt, nt)
    @test rel_err(sol.q, refx) < 1E-7


    # Test integration with exact Jacobian

    int = IntegratorFIRK(ode, TableauGLRK(1), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 5E-4

    int = IntegratorFIRK(ode, TableauGLRK(2), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 1E-7

    int = IntegratorFIRK(ode, TableauGLRK(3), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorFIRK(ode, TableauGLRK(4), Δt; exact_jacobian=true)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, refx) < 1E-15

end


@testset "$(rpad("Partitioned Runge-Kutta integrators",80))" begin

    psol = integrate(pode, TableauSymplecticEulerA(), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-2
    @test rel_err(psol.p, refp) < 1E-3

    psol = integrate(pode, TableauSymplecticEulerB(), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-2
    @test rel_err(psol.p, refp) < 1E-3

    psol = integrate(pode, TableauEPRK(:prk4, 4, TableauRK4().q), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-7
    @test rel_err(psol.p, refp) < 5E-7

    psol = integrate(pode, TableauIPRK(:pglrk, 2, CoefficientsGLRK(1)), Δt, nt)
    @test rel_err(psol.q, refq) < 5E-4
    @test rel_err(psol.p, refp) < 5E-4

    psol = integrate(pode, TableauIPRK(:pglrk, 4, CoefficientsGLRK(2)), Δt, nt)
    @test rel_err(psol.q, refq) < 1E-7
    @test rel_err(psol.p, refp) < 1E-7

    psol = integrate(pode, TableauIPRK(:pglrk, 6, CoefficientsGLRK(3)), Δt, nt)
    @test rel_err(psol.q, refq) < 1E-11
    @test rel_err(psol.p, refp) < 1E-11

end


@testset "$(rpad("Special Runge-Kutta integrators",80))" begin

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(2), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 4E-8

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(3), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 2E-12

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(4), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 8E-16

end


using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.Oscillator
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.Oscillator: Δt, nt, refx, refq, refp

ode  = oscillator_ode()
pode = oscillator_pode()
sode = oscillator_sode()
iode = oscillator_iode()
idae = oscillator_idae()
pdae = oscillator_pdae()


@testset "$(rpad("Runge-Kutta integrators",80))" begin

    @test typeof(Integrator(ode, getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
    @test typeof(Integrator(ode, getTableauCrouzeix(), Δt)) <: IntegratorDIRK
    @test typeof(Integrator(ode, getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK
    @test typeof(Integrator(ode, getTableauGLRK(1), Δt)) <: IntegratorFIRK

end


@testset "$(rpad("Explicit Runge-Kutta integrators",80))" begin

    int = Integrator(ode, getTableauExplicitEuler(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-2

    int = Integrator(ode, getTableauExplicitMidpoint(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 1E-3

    int = Integrator(ode, getTableauERK4(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-7

end


@testset "$(rpad("Implicit Runge-Kutta integrators",80))" begin

    int = Integrator(ode, getTableauCrouzeix(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-5

    int = Integrator(ode, getTableauImplicitEuler(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-2

    int = Integrator(ode, getTableauImplicitMidpoint(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-4

    int = Integrator(ode, getTableauGLRK(1), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-4

    int = Integrator(ode, getTableauGLRK(2), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 1E-7

    int = Integrator(ode, getTableauGLRK(3), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(ode, getTableauGLRK(4), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 1E-15

    int = Integrator(ode, getTableauGLRK(5), Δt)

    @test rel_err(sol.q, refx) < 1E-15

    int = Integrator(ode, getTableauGLRK(6), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 1E-15

    int = Integrator(ode, getTableauGLRK(7), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 5E-15

    int = Integrator(ode, getTableauSRK3(), Δt)
    sol = integrate(int, nt)

    @test rel_err(sol.q, refx) < 1E-7

end


@testset "$(rpad("Partitioned Runge-Kutta integrators",80))" begin

    pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
    psol = integrate(pint, nt)

    @test rel_err(psol.q, refq) < 5E-2
    @test rel_err(psol.p, refp) < 1E-3

    pint = Integrator(pode, getTableauSymplecticEulerB(), Δt)
    psol = integrate(pint, nt)

    @test rel_err(psol.q, refq) < 5E-2
    @test rel_err(psol.p, refp) < 1E-3

    pint = Integrator(pode, TableauEPRK(:prk4, 4, getTableauERK4().q), Δt)
    psol = integrate(pint, nt)

    @test rel_err(psol.q, refq) < 5E-7
    @test rel_err(psol.p, refp) < 5E-7

    pint = Integrator(pode, TableauIPRK(:pglrk, 2, getCoefficientsGLRK(1)), Δt)
    psol = integrate(pint, nt)

    @test rel_err(psol.q, refq) < 5E-4
    @test rel_err(psol.p, refp) < 5E-4

    pint = Integrator(pode, TableauIPRK(:pglrk, 4, getCoefficientsGLRK(2)), Δt)
    psol = integrate(pint, nt)

    @test rel_err(psol.q, refq) < 1E-7
    @test rel_err(psol.p, refp) < 1E-7

    pint = Integrator(pode, TableauIPRK(:pglrk, 6, getCoefficientsGLRK(3)), Δt)
    psol = integrate(pint, nt)

    @test rel_err(psol.q, refq) < 1E-11
    @test rel_err(psol.p, refp) < 1E-11

end


@testset "$(rpad("Special Runge-Kutta integrators",80))" begin

    pgint = IntegratorPGLRK(iode, getCoefficientsPGLRK(2), Δt)
    pgsol = integrate(pgint, nt)

    @test rel_err(pgsol.q, refx) < 1E-5

end

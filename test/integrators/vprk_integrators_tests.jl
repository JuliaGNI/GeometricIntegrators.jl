
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: Δt, nt

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

ode  = lotka_volterra_2d_ode()
iode = lotka_volterra_2d_iode()
vode = lotka_volterra_2d_vode()

int  = IntegratorFIRK(ode, getTableauGLRK(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[:,end]


@testset "$(rpad("VPRK integrators",80))" begin

    sol = integrate(iode, TableauVPRK(:pglrk, 2, getCoefficientsGLRK(1), -1), Δt, nt)
    @test rel_err(sol.q, refx) < 2E-6

    sol = integrate(iode, TableauVPRK(:pglrk, 4, getCoefficientsGLRK(2), +1), Δt, nt)
    @test rel_err(sol.q, refx) < 8E-7

    sol = integrate(iode, TableauVPRK(:pglrk, 6, getCoefficientsGLRK(3), -1), Δt, nt)
    @test rel_err(sol.q, refx) < 4E-12

    sol = integrate(iode, getTableauVPLobIIIA2(), Δt, nt)
    @test rel_err(sol.q, refx) < 4E-6

    sol = integrate(iode, getTableauVPLobIIIA3(), Δt, nt)
    @test rel_err(sol.q, refx) < 8E-7

    sol = integrate(iode, getTableauVPLobIIIA4(), Δt, nt)
    @test rel_err(sol.q, refx) < 3E-11

    sol = integrate(iode, getTableauVPLobIIIB2(), Δt, nt)
    @test rel_err(sol.q, refx) < 2E-6

    sol = integrate(iode, getTableauVPLobIIIB3(), Δt, nt)
    @test rel_err(sol.q, refx) < 8E-7

    sol = integrate(iode, getTableauVPLobIIIB4(), Δt, nt)
    @test rel_err(sol.q, refx) < 2E-11

end


@testset "$(rpad("VPRK integrators with standard projection",80))" begin

    int = IntegratorVPRKpStandard(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVPRKpStandard(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVPRKpStandard(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 5E-16

end


@testset "$(rpad("VPRK integrators with symplectic projection",80))" begin

    int = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 8E-12

end


@testset "$(rpad("VPRK integrators with symmetric projection",80))" begin

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

end


@testset "$(rpad("VPRK integrators with midpoint projection",80))" begin

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

end


@testset "$(rpad("VPRK integrators with internal projection",80))" begin

    int = IntegratorVPRKpInternal(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = IntegratorVPRKpInternal(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVPRKpInternal(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-12

    int = IntegratorVPRKpInternal(iode, getTableauVPGLRK(4), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 8E-16

end


@testset "$(rpad("VPRK integrators with projection on secondary constraint",80))" begin

    int = IntegratorVPRKpSecondary(vode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = IntegratorVPRKpSecondary(vode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 8E-7

    int = IntegratorVPRKpSecondary(vode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-12

end


@testset "$(rpad("VPRK integrators with variational projection",80))" begin

    intV1 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(1), Δt)
    solV1 = integrate(iode, intV1, nt)
    @test rel_err(solV1.q, refx) < 8E-7

    intV2 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(2), Δt)
    solV2 = integrate(iode, intV2, nt)
    @test rel_err(solV2.q, refx) < 8E-8

    intV3 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(3), Δt)
    solV3 = integrate(iode, intV3, nt)
    @test rel_err(solV3.q, refx) < 1E-11

    intQ1 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(1), Δt)
    solQ1 = integrate(iode, intQ1, nt)
    @test rel_err(solQ1.q, refx) < 4E-5

    intQ2 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(2), Δt)
    solQ2 = integrate(iode, intQ2, nt)
    @test rel_err(solQ2.q, refx) < 2E-4

    intQ3 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(3), Δt)
    solQ3 = integrate(iode, intQ3, nt)
    @test rel_err(solQ3.q, refx) < 1E-8

    intP1 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(1), Δt)
    solP1 = integrate(iode, intP1, nt)
    @test rel_err(solP1.q, refx) < 8E-7

    intP2 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(2), Δt)
    solP2 = integrate(iode, intP2, nt)
    @test rel_err(solP2.q, refx) < 8E-8

    intP3 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(3), Δt)
    solP3 = integrate(iode, intP3, nt)
    @test rel_err(solP3.q, refx) < 1E-11

    @test rel_err(solV1.q, solP1.q[:,end]) == 0
    @test rel_err(solV2.q, solP2.q[:,end]) == 0
    @test rel_err(solV3.q, solP3.q[:,end]) == 0

end


@testset "$(rpad("VPRK integrators with projection on Runge-Kutta tableau",80))" begin

    int = IntegratorVPRKpTableau(iode, getCoefficientsPGLRK(5), Δt*5)
    sol = integrate(iode, int, div(nt,5))
    @test rel_err(sol.q, refx) < 8E-12

    int = IntegratorVPRKpTableau(iode, getCoefficientsPGLRK(6), Δt*5)
    sol = integrate(iode, int, div(nt,5))
    @test rel_err(sol.q, refx) < 4E-14

end


@testset "$(rpad("Degenerate symplectic partitioned Runge-Kutta methods",80))" begin

    int = IntegratorVPRKdegenerate(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 2E-5

    int = IntegratorVPRKdegenerate(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-7

    int = IntegratorVPRKdegenerate(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 2E-10

end


@testset "$(rpad("VSPRK integrators with Legendre projection",80))" begin

    int = IntegratorVPRKpLegendre(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVPRKpLegendre(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVPRKpLegendre(iode, getTableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

end

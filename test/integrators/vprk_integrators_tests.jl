
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: Δt, nt

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

iode = lotka_volterra_2d_iode()
vode = lotka_volterra_2d_vode()

int  = IntegratorFIRK(lotka_volterra_2d_ode(), getTableauGLRK(8), Δt)
sol  = integrate(int, nt)
refx = sol.q[:,end]


@testset "$(rpad("VPRK integrators",80))" begin

    iint = Integrator(iode, TableauVPRK(:pglrk, 2, getCoefficientsGLRK(1), -1), Δt)
    isol = integrate(iint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 2E-6

    iint = Integrator(iode, TableauVPRK(:pglrk, 4, getCoefficientsGLRK(2), +1), Δt)
    isol = integrate(iint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 8E-7

    iint = Integrator(iode, TableauVPRK(:pglrk, 6, getCoefficientsGLRK(3), -1), Δt)
    isol = integrate(iint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-12

    vint = Integrator(iode, getTableauVPLobIIIA2(), Δt)
    vsol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(vsol.q, refx) < 4E-6

    vint = Integrator(iode, getTableauVPLobIIIA3(), Δt)
    vsol = integrate(vint, nt)

    # println(rel_err(vsol.q, refx))
    @test rel_err(vsol.q, refx) < 8E-7

    vint = Integrator(iode, getTableauVPLobIIIA4(), Δt)
    vsol = integrate(vint, nt)

    # println(rel_err(vsol.q, refx))
    @test rel_err(vsol.q, refx) < 3E-11

    vint = Integrator(iode, getTableauVPLobIIIB2(), Δt)
    vsol = integrate(vint, nt)

    # println(rel_err(vsol.q, refx))
    @test rel_err(vsol.q, refx) < 2E-6

    vint = Integrator(iode, getTableauVPLobIIIB3(), Δt)
    vsol = integrate(vint, nt)

    # println(rel_err(vsol.q, refx))
    @test rel_err(vsol.q, refx) < 8E-7

    vint = Integrator(iode, getTableauVPLobIIIB4(), Δt)
    vsol = integrate(vint, nt)

    # println(rel_err(vsol.q, refx))
    @test rel_err(vsol.q, refx) < 2E-11

end


@testset "$(rpad("VPRK integrators with standard projection",80))" begin

    vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(1), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-6

    vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(2), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-11

    vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(3), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 5E-16

end


@testset "$(rpad("VPRK integrators with symplectic projection",80))" begin

    vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(1), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-6

    vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(2), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-11

    vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(3), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 8E-12

end


@testset "$(rpad("VPRK integrators with symmetric projection",80))" begin

    vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-6

    vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(2), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-11

    vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(3), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-16

end


@testset "$(rpad("VPRK integrators with midpoint projection",80))" begin

    vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(1), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-6

    vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(2), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-11

    vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(3), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-16

end


@testset "$(rpad("VPRK integrators with internal projection",80))" begin

    vint = IntegratorVPRKpInternal(iode, getTableauVPGLRK(1), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 2E-6

    vint = IntegratorVPRKpInternal(iode, getTableauVPGLRK(2), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 1E-11

    vint = IntegratorVPRKpInternal(iode, getTableauVPGLRK(3), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-12

    vint = IntegratorVPRKpInternal(iode, getTableauVPGLRK(4), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 8E-16

end


@testset "$(rpad("VPRK integrators with projection on secondary constraint",80))" begin

    vint = IntegratorVPRKpSecondary(vode, getTableauVPGLRK(1), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 2E-6

    vint = IntegratorVPRKpSecondary(vode, getTableauVPGLRK(2), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 8E-7

    vint = IntegratorVPRKpSecondary(vode, getTableauVPGLRK(3), Δt)
    isol = integrate(vint, nt)

    # println(rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-12

end


@testset "$(rpad("VPRK integrators with variational projection",80))" begin

    vintV1 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(1), Δt)
    isolV1 = integrate(vintV1, nt)

    # println(rel_err(isolV1.q, refx))
    @test rel_err(isolV1.q, refx) < 8E-7

    vintV2 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(2), Δt)
    isolV2 = integrate(vintV2, nt)

    # println(rel_err(isolV2.q, refx))
    @test rel_err(isolV2.q, refx) < 8E-8

    vintV3 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(3), Δt)
    isolV3 = integrate(vintV3, nt)

    # println(rel_err(isolV3.q, refx))
    @test rel_err(isolV3.q, refx) < 1E-11

    vintQ1 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(1), Δt)
    isolQ1 = integrate(vintQ1, nt)

    # println(rel_err(isolQ1.q, refx))
    @test rel_err(isolQ1.q, refx) < 4E-5

    vintQ2 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(2), Δt)
    isolQ2 = integrate(vintQ2, nt)

    # println(rel_err(isolQ2.q, refx))
    @test rel_err(isolQ2.q, refx) < 2E-4

    vintQ3 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(3), Δt)
    isolQ3 = integrate(vintQ3, nt)

    # println(rel_err(isolQ3.q, refx))
    @test rel_err(isolQ3.q, refx) < 1E-8

    vintP1 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(1), Δt)
    isolP1 = integrate(vintP1, nt)

    # println(rel_err(isolP1.q, refx))
    @test rel_err(isolP1.q, refx) < 8E-7

    vintP2 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(2), Δt)
    isolP2 = integrate(vintP2, nt)

    # println(rel_err(isolP2.q, refx))
    @test rel_err(isolP2.q, refx) < 8E-8

    vintP3 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(3), Δt)
    isolP3 = integrate(vintP3, nt)

    # println(rel_err(isolP3.q, refx))
    @test rel_err(isolP3.q, refx) < 1E-11

    @test rel_err(isolV1.q, isolP1.q[:,end]) == 0
    @test rel_err(isolV2.q, isolP2.q[:,end]) == 0
    @test rel_err(isolV3.q, isolP3.q[:,end]) == 0

end


@testset "$(rpad("VPRK integrators with projection on Runge-Kutta tableau",80))" begin

    vint = IntegratorVPRKpTableau(iode, getCoefficientsPGLRK(5), Δt*5)
    isol = integrate(vint, div(nt,5))

    # println("error = ", rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 8E-12

    vint = IntegratorVPRKpTableau(iode, getCoefficientsPGLRK(6), Δt*5)
    isol = integrate(vint, div(nt,5))

    # println("error = ", rel_err(isol.q, refx))
    @test rel_err(isol.q, refx) < 4E-14

end

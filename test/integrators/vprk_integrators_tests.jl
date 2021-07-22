
using GeometricBase.Config
using GeometricBase.Utils
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Tableaus
using GeometricProblems.LotkaVolterra2d
using SimpleSolvers
using Test

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 10
const q₀ = [1.0, 1.0]
const parameters = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

ode  = lotka_volterra_2d_ode(q₀; params=parameters)
iode = lotka_volterra_2d_iode(q₀; params=parameters)
lode = lotka_volterra_2d_lode(q₀; params=parameters)
ldae = lotka_volterra_2d_ldae(q₀; params=parameters)

int  = IntegratorFIRK(ode, TableauGauss(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[end]


@testset "$(rpad("VPRK integrators",80))" begin

    sol = integrate(iode, TableauVPRK(:pglrk, 2, TableauGauss(1), -1), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 2E-6

    sol = integrate(iode, TableauVPRK(:pglrk, 4, TableauGauss(2), +1), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 8E-7

    sol = integrate(iode, TableauVPRK(:pglrk, 6, TableauGauss(3), -1), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-12

    sol = integrate(iode, TableauVPLobattoIIIA(2), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-6

    sol = integrate(iode, TableauVPLobattoIIIA(3), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 8E-7

    sol = integrate(iode, TableauVPLobattoIIIA(4), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 3E-11

    sol = integrate(iode, TableauVPLobattoIIIB(2), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 2E-6

    sol = integrate(iode, TableauVPLobattoIIIB(3), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 8E-7

    sol = integrate(iode, TableauVPLobattoIIIB(4), Δt, nt)
    @test relative_maximum_error(sol.q, refx) < 2E-11

end


@testset "$(rpad("VPRK integrators with standard projection",80))" begin

    int = IntegratorVPRKpStandard(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-6

    int = IntegratorVPRKpStandard(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorVPRKpStandard(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-15

end


@testset "$(rpad("VPRK integrators with symplectic projection",80))" begin

    int = IntegratorVPRKpSymplectic(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-6

    int = IntegratorVPRKpSymplectic(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorVPRKpSymplectic(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 8E-12

end


@testset "$(rpad("VPRK integrators with symmetric projection",80))" begin

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-6

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-16

end


@testset "$(rpad("VPRK integrators with midpoint projection",80))" begin

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-6

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-16

end


@testset "$(rpad("VPRK integrators with internal projection",80))" begin

    int = IntegratorVPRKpInternal(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 2E-6

    int = IntegratorVPRKpInternal(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorVPRKpInternal(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-12

    int = IntegratorVPRKpInternal(iode, TableauVPGLRK(4), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 8E-16

end


@testset "$(rpad("VPRK integrators with projection on secondary constraint",80))" begin

    # TODO: reactivate

    # int = IntegratorVPRKpSecondary(ldae, TableauVPGLRK(1), Δt)
    # sol = integrate(ldae, int, nt)
    # @test relative_maximum_error(sol.q, refx) < 2E-6

    # int = IntegratorVPRKpSecondary(ldae, TableauVPGLRK(2), Δt)
    # sol = integrate(ldae, int, nt)
    # @test relative_maximum_error(sol.q, refx) < 8E-7

    # int = IntegratorVPRKpSecondary(ldae, TableauVPGLRK(3), Δt)
    # sol = integrate(ldae, int, nt)
    # @test relative_maximum_error(sol.q, refx) < 4E-12

end


@testset "$(rpad("VPRK integrators with variational projection",80))" begin

    intV1 = IntegratorVPRKpVariational(iode, TableauVPGLRK(1), Δt)
    solV1 = integrate(iode, intV1, nt)
    @test relative_maximum_error(solV1.q, refx) < 8E-7

    intV2 = IntegratorVPRKpVariational(iode, TableauVPGLRK(2), Δt)
    solV2 = integrate(iode, intV2, nt)
    @test relative_maximum_error(solV2.q, refx) < 8E-8

    intV3 = IntegratorVPRKpVariational(iode, TableauVPGLRK(3), Δt)
    solV3 = integrate(iode, intV3, nt)
    @test relative_maximum_error(solV3.q, refx) < 1E-11

    intQ1 = IntegratorVPRKpVariationalQ(iode, TableauVPGLRK(1), Δt)
    solQ1 = integrate(iode, intQ1, nt)
    @test relative_maximum_error(solQ1.q, refx) < 4E-5

    intQ2 = IntegratorVPRKpVariationalQ(iode, TableauVPGLRK(2), Δt)
    solQ2 = integrate(iode, intQ2, nt)
    @test relative_maximum_error(solQ2.q, refx) < 2E-4

    intQ3 = IntegratorVPRKpVariationalQ(iode, TableauVPGLRK(3), Δt)
    solQ3 = integrate(iode, intQ3, nt)
    @test relative_maximum_error(solQ3.q, refx) < 1E-8

    intP1 = IntegratorVPRKpVariationalP(iode, TableauVPGLRK(1), Δt)
    solP1 = integrate(iode, intP1, nt)
    @test relative_maximum_error(solP1.q, refx) < 8E-7

    intP2 = IntegratorVPRKpVariationalP(iode, TableauVPGLRK(2), Δt)
    solP2 = integrate(iode, intP2, nt)
    @test relative_maximum_error(solP2.q, refx) < 8E-8

    intP3 = IntegratorVPRKpVariationalP(iode, TableauVPGLRK(3), Δt)
    solP3 = integrate(iode, intP3, nt)
    @test relative_maximum_error(solP3.q, refx) < 1E-11

    @test relative_maximum_error(solV1.q, solP1.q[end]) == 0
    @test relative_maximum_error(solV2.q, solP2.q[end]) == 0
    @test relative_maximum_error(solV3.q, solP3.q[end]) == 0

end


@testset "$(rpad("VPRK integrators with projection on Runge-Kutta tableau",80))" begin

    int = IntegratorVPRKpTableau(iode, CoefficientsPGLRK(5), Δt*5)
    sol = integrate(iode, int, div(nt,5))
    @test relative_maximum_error(sol.q, refx) < 8E-12

    int = IntegratorVPRKpTableau(iode, CoefficientsPGLRK(6), Δt*5)
    sol = integrate(iode, int, div(nt,5))
    @test relative_maximum_error(sol.q, refx) < 4E-14

end


@testset "$(rpad("Degenerate symplectic partitioned Runge-Kutta methods",80))" begin

    int = IntegratorVPRKdegenerate(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 2E-5

    int = IntegratorVPRKdegenerate(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 4E-7

    int = IntegratorVPRKdegenerate(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 2E-10

end


@testset "$(rpad("VSPRK integrators with Legendre projection",80))" begin

    int = IntegratorVPRKpLegendre(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-6

    int = IntegratorVPRKpLegendre(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 1E-11

    int = IntegratorVPRKpLegendre(iode, TableauVPGLRK(3), Δt)
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, refx) < 8E-16

end

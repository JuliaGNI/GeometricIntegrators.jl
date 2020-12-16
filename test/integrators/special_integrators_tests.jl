
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem
using GeometricIntegrators.Utils
using SimpleSolvers
using Test

using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: Δt, nt

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

ode  = lotka_volterra_2d_ode()
vode = lotka_volterra_2d_vode()
int  = IntegratorFIRK(ode, TableauGLRK(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[:,end]


@testset "$(rpad("Special integrators",80))" begin

    flint = IntegratorFLRK(vode, TableauGLRK(2), Δt)
    flsol = integrate(vode, flint, nt)
    @test rel_err(flsol.q, refx) < 4E-12

    flint = IntegratorFLRK(vode, TableauGLRK(3), Δt)
    flsol = integrate(vode, flint, nt)
    @test rel_err(flsol.q, refx) < 2E-16

    flint = IntegratorFLRK(vode, TableauGLRK(4), Δt)
    flsol = integrate(vode, flint, nt)
    @test rel_err(flsol.q, refx) < 2E-16


    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(2), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 4E-12

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(3), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 2E-16

    pgint = IntegratorPGLRK(ode, CoefficientsPGLRK(4), Δt)
    pgsol = integrate(ode, pgint, nt)
    @test rel_err(pgsol.q, refx) < 2E-16

end

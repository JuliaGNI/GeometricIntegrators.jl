
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricIntegrators.Utils
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
vode = lotka_volterra_2d_vode(q₀; params=parameters)
int  = IntegratorFIRK(ode, TableauGauss(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[end]


@testset "$(rpad("Special integrators",80))" begin

    flint = IntegratorFLRK(vode, TableauGauss(2), Δt)
    flsol = integrate(vode, flint, nt)
    @test rel_err(flsol.q, refx) < 4E-12

    flint = IntegratorFLRK(vode, TableauGauss(3), Δt)
    flsol = integrate(vode, flint, nt)
    @test rel_err(flsol.q, refx) < 2E-16

    flint = IntegratorFLRK(vode, TableauGauss(4), Δt)
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

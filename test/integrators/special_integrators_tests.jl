
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

ode = lotka_volterra_2d_ode()

int  = IntegratorFIRK(ode, getTableauGLRK(8), Δt)
sol  = integrate(int, nt)
refx = sol.q[:,end]


@testset "$(rpad("Special integrators",80))" begin

    pgint = IntegratorPGLRK(ode, getCoefficientsPGLRK(2), Δt)
    pgsol = integrate(pgint, nt)

    # println(rel_err(pgsol.q, refx))
    @test rel_err(pgsol.q, refx) < 4E-12

    pgint = IntegratorPGLRK(ode, getCoefficientsPGLRK(3), Δt)
    pgsol = integrate(pgint, nt)

    # println(rel_err(pgsol.q, refx))
    @test rel_err(pgsol.q, refx) < 2E-16

    pgint = IntegratorPGLRK(ode, getCoefficientsPGLRK(4), Δt)
    pgsol = integrate(pgint, nt)

    # println(rel_err(pgsol.q, refx))
    @test rel_err(pgsol.q, refx) < 2E-16

end

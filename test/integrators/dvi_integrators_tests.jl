
using GeometricBase.Config
using GeometricBase.Utils
using GeometricIntegrators.Integrators
using GeometricIntegrators.Solutions
using GeometricIntegrators.Tableaus
using GeometricProblems.LotkaVolterra2d
using SimpleSolvers
using Test

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 10
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

ode  = lotka_volterra_2d_ode(q₀; parameters=params, tstep=Δt)
iode = lotka_volterra_2d_iode(q₀; parameters=params, tstep=Δt)
lode = lotka_volterra_2d_lode(q₀; parameters=params, tstep=Δt)

int  = IntegratorFIRK(ode, TableauGauss(8))
sol  = integrate(ode, int, nt)
refx = sol.q[end]


@testset "$(rpad("1st Order DVIs",80))" begin

    sol = Solution(iode, nt)
    int = IntegratorDVIA(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 1E-1

    sol = Solution(iode, nt)
    int = IntegratorDVIB(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 1E-1

end


@testset "$(rpad("2nd Order Centred DVIs",80))" begin

    sol = Solution(iode, nt)
    int = IntegratorCMDVI(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 4E-3

    sol = Solution(iode, nt)
    int = IntegratorCTDVI(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 4E-3

end

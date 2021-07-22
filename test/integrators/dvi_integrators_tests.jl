
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
const parameters = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

ode  = lotka_volterra_2d_ode(q₀; params=parameters)
iode = lotka_volterra_2d_iode(q₀; params=parameters)
lode = lotka_volterra_2d_lode(q₀; params=parameters)

int  = IntegratorFIRK(ode, TableauGauss(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[end]


@testset "$(rpad("1st Order DVIs",80))" begin

    sol = Solution(iode, Δt, nt)
    int = IntegratorDVIA(iode, Δt)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 1E-1

    sol = Solution(iode, Δt, nt)
    int = IntegratorDVIB(iode, Δt)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 1E-1

end


@testset "$(rpad("2nd Order Centred DVIs",80))" begin

    sol = Solution(iode, Δt, nt)
    int = IntegratorCMDVI(iode, Δt)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 4E-3

    sol = Solution(iode, Δt, nt)
    int = IntegratorCTDVI(iode, Δt)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 4E-3

end

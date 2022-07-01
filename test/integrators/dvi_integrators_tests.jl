using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using SimpleSolvers
using Test

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt*nt)

ode  = lotka_volterra_2d_ode(q₀; tspan=tspan, tstep=Δt, parameters=params)
iode = lotka_volterra_2d_iode(q₀; tspan=tspan, tstep=Δt, parameters=params)
lode = lotka_volterra_2d_lode(q₀; tspan=tspan, tstep=Δt, parameters=params)

int  = IntegratorFIRK(ode, TableauGauss(8))
sol  = integrate(ode, int)
refx = sol.q[end]


@testset "$(rpad("1st Order DVIs",80))" begin

    sol = Solution(iode)
    int = IntegratorDVIA(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 1E-1

    sol = Solution(iode)
    int = IntegratorDVIB(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 1E-1

end


@testset "$(rpad("2nd Order Centred DVIs",80))" begin

    sol = Solution(iode)
    int = IntegratorCMDVI(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 4E-3

    sol = Solution(iode)
    int = IntegratorCTDVI(iode)
    integrate!(int, sol)
    @test relative_maximum_error(sol.q, refx) < 4E-3

end

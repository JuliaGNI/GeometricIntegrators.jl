using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using SimpleSolvers
using Test


const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt*nt)

ode  = lotka_volterra_2d_ode(q₀; tspan=tspan, tstep=Δt, parameters=params)
iode = lotka_volterra_2d_iode(q₀; tspan=tspan, tstep=Δt, parameters=params)
lode = lotka_volterra_2d_lode(q₀; tspan=tspan, tstep=Δt, parameters=params)

sol  = integrate(ode, Gauss(8))
reference_solution = sol.q[end]


@testset "$(rpad("1st Order DVIs",80))" begin

    sol = integrate(lode, DVIA())
    @test relative_maximum_error(sol.q, reference_solution) < 1E-1

    sol = integrate(lode, DVIB())
    @test relative_maximum_error(sol.q, reference_solution) < 1E-1

end


@testset "$(rpad("2nd Order Centred DVIs",80))" begin

    sol = integrate(lode, CMDVI())
    @test relative_maximum_error(sol.q, reference_solution) < 4E-3

    sol = integrate(lode, CTDVI())
    @test relative_maximum_error(sol.q, reference_solution) < 4E-3

end


@testset "$(rpad("Degenerate Variational Runge-Kutta integrators",80))" begin

    sol = integrate(lode, DVRK(Gauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-5

    sol = integrate(lode, DVRK(Gauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-7

    sol = integrate(lode, DVRK(Gauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-10

    sol = integrate(lode, DVRK(Gauss(4)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-13

end

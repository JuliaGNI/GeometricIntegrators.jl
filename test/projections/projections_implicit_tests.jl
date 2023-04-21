using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using Test


const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt*nt)

ode  = lotka_volterra_2d_ode(q₀; tspan=tspan, tstep=Δt, parameters=params)
iode = lotka_volterra_2d_iode(q₀; tspan=tspan, tstep=Δt, parameters=params)

sol  = integrate(ode, Gauss(8))

reference_solution = sol.q[end]


@testset "$(rpad("Post-projection with Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, PostProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-4

    sol = integrate(iode, PostProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-8

    sol = integrate(iode, PostProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-12

    sol = integrate(iode, PostProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-15

end


@testset "$(rpad("Midpoint Projection with Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, MidpointProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-4

    sol = integrate(iode, MidpointProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-8

    sol = integrate(iode, MidpointProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-12

    sol = integrate(iode, MidpointProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-15

end

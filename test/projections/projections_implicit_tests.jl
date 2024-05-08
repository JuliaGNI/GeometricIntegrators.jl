using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using Test


const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt*nt)

ode  = odeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
iode = iodeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
ref  = integrate(ode, Gauss(8))


@testset "$(rpad("Post-projection with Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, PostProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(iode, PostProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-8

    sol = integrate(iode, PostProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-12

    sol = integrate(iode, PostProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15

end


@testset "$(rpad("Midpoint Projection with Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, MidpointProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(iode, MidpointProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-8

    sol = integrate(iode, MidpointProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-12

    sol = integrate(iode, MidpointProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15

end


@testset "$(rpad("Symmetric Projection with Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, SymmetricProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(iode, SymmetricProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-8

    sol = integrate(iode, SymmetricProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-12

    sol = integrate(iode, SymmetricProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15

end

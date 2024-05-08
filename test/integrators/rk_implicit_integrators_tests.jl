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
lode = lodeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)

ref  = integrate(ode, Gauss(8))


@testset "$(rpad("Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, Gauss(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-5

    sol = integrate(iode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-7

    sol = integrate(iode, Gauss(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-10

    sol = integrate(iode, Gauss(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-13


    sol = integrate(iode, IRK(Gauss(1); implicit_update = true))
    @test relative_maximum_error(sol.q, ref.q) < 4E-5

    sol = integrate(iode, IRK(Gauss(2); implicit_update = true))
    @test relative_maximum_error(sol.q, ref.q) < 8E-7

    sol = integrate(iode, IRK(Gauss(3); implicit_update = true))
    @test relative_maximum_error(sol.q, ref.q) < 4E-10

    sol = integrate(iode, IRK(Gauss(4); implicit_update = true))
    @test relative_maximum_error(sol.q, ref.q) < 2E-13

end


@testset "$(rpad("Partitioned Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, PartitionedGauss(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-5

    sol = integrate(iode, PartitionedGauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-7

    sol = integrate(iode, PartitionedGauss(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-10

    sol = integrate(iode, PartitionedGauss(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-13

end


# @testset "$(rpad("Formal Lagrangian Runge-Kutta integrators",80))" begin

#     flint = IntegratorFLRK(lode, TableauGauss(2))
#     flsol = integrate(lode, flint)
#     @test relative_maximum_error(flsol.q, ref.q) < 4E-12

#     flint = IntegratorFLRK(lode, TableauGauss(3))
#     flsol = integrate(lode, flint)
#     @test relative_maximum_error(flsol.q, ref.q) < 4E-16

#     flint = IntegratorFLRK(lode, TableauGauss(4))
#     flsol = integrate(lode, flint)
#     @test relative_maximum_error(flsol.q, ref.q) < 8E-16

# end

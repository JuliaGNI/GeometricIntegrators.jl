using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using Test


const Δt = 0.01
const nt = 10
const q₀ = [1.0, 1.0]
const tspan = (0.0, Δt * nt)
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

ode = lotka_volterra_2d_ode(q₀; tspan=tspan, tstep=Δt, parameters=params)
iode = lotka_volterra_2d_iode(q₀; tspan=tspan, tstep=Δt, parameters=params)
lode = lotka_volterra_2d_lode(q₀; tspan=tspan, tstep=Δt, parameters=params)
ldae = lotka_volterra_2d_ldae(q₀; tspan=tspan, tstep=Δt, parameters=params)

sol = integrate(ode, Gauss(8))
reference_solution = sol.q[end]


@testset "$(rpad("VPRK integrators without projection",80))" begin

    sol = integrate(iode, VPRKGauss(1))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-6

    sol = integrate(iode, VPRKGauss(2))
    @test relative_maximum_error(sol.q, reference_solution) < 8E-7

    sol = integrate(iode, VPRKGauss(3))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-12

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(2))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-6

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(3))
    @test relative_maximum_error(sol.q, reference_solution) < 8E-7

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(4))
    @test relative_maximum_error(sol.q, reference_solution) < 3E-11

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(2))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-6

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(3))
    @test relative_maximum_error(sol.q, reference_solution) < 8E-7

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(4))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-11

end


@testset "$(rpad("VPRK integrators with standard projection",80))" begin

    sol = integrate(iode, PostProjection(VPRKGauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-6

    sol = integrate(iode, PostProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-11

    sol = integrate(iode, PostProjection(VPRKGauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-15

end


@testset "$(rpad("VPRK integrators with symplectic projection",80))" begin

    sol = integrate(iode, SymplecticProjection(VPRKGauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-6

    sol = integrate(iode, SymplecticProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-11

    sol = integrate(iode, SymplecticProjection(VPRKGauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-15

end


@testset "$(rpad("VPRK integrators with midpoint projection",80))" begin

    sol = integrate(iode, MidpointProjection(VPRKGauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-6

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-11

    sol = integrate(iode, MidpointProjection(VPRKGauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-15

end


@testset "$(rpad("VPRK integrators with symmetric projection",80))" begin

    sol = integrate(iode, SymmetricProjection(VPRKGauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-6

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 1E-11

    sol = integrate(iode, SymmetricProjection(VPRKGauss(3)))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-15

end


# @testset "$(rpad("VPRK integrators with internal projection",80))" begin

#     sol = integrate(iode, VPRKpInternal(VPRKGauss(1)))
#     @test relative_maximum_error(sol.q, reference_solution) < 2E-6

#     sol = integrate(iode, VPRKpInternal(VPRKGauss(2)))
#     @test relative_maximum_error(sol.q, reference_solution) < 1E-11

#     sol = integrate(iode, VPRKpInternal(VPRKGauss(3)))
#     @test relative_maximum_error(sol.q, reference_solution) < 4E-12

#     sol = integrate(iode, VPRKpInternal(VPRKGauss(4)))
#     @test relative_maximum_error(sol.q, reference_solution) < 4E-15

# end


# @testset "$(rpad("VPRK integrators with projection on secondary constraint",80))" begin

#     # TODO: reactivate

#     # sol = integrate(ldae, VPRKpSecondary(VPRKGauss(1)))
#     # @test relative_maximum_error(sol.q, reference_solution) < 2E-6

#     # sol = integrate(ldae, VPRKpSecondary(VPRKGauss(2)))
#     # @test relative_maximum_error(sol.q, reference_solution) < 8E-7

#     # sol = integrate(ldae, VPRKpSecondary(VPRKGauss(3)))
#     # @test relative_maximum_error(sol.q, reference_solution) < 4E-12

# end


# @testset "$(rpad("VPRK integrators with variational projection",80))" begin

#     solV1 = integrate(iode, VPRKpVariational(VPRKGauss(1)))
#     @test relative_maximum_error(solV1.q, reference_solution) < 8E-7

#     solV2 = integrate(iode, VPRKpVariational(VPRKGauss(2)))
#     @test relative_maximum_error(solV2.q, reference_solution) < 8E-8

#     solV3 = integrate(iode, VPRKpVariational(VPRKGauss(3)))
#     @test relative_maximum_error(solV3.q, reference_solution) < 1E-11

#     solQ1 = integrate(iode, VPRKpVariationalQ(VPRKGauss(1)))
#     @test relative_maximum_error(solQ1.q, reference_solution) < 4E-5

#     solQ2 = integrate(iode, VPRKpVariationalQ(VPRKGauss(2)))
#     @test relative_maximum_error(solQ2.q, reference_solution) < 2E-4

#     solQ3 = integrate(iode, VPRKpVariationalQ(VPRKGauss(3)))
#     @test relative_maximum_error(solQ3.q, reference_solution) < 1E-8

#     solP1 = integrate(iode, VPRKpVariationalP(VPRKGauss(1)))
#     @test relative_maximum_error(solP1.q, reference_solution) < 8E-7

#     solP2 = integrate(iode, VPRKpVariationalP(VPRKGauss(2)))
#     @test relative_maximum_error(solP2.q, reference_solution) < 8E-8

#     solP3 = integrate(iode, VPRKpVariationalP(VPRKGauss(3)))
#     @test relative_maximum_error(solP3.q, reference_solution) < 1E-11

#     @test relative_maximum_error(solV1.q, solP1.q[end]) == 0
#     @test relative_maximum_error(solV2.q, solP2.q[end]) == 0
#     @test relative_maximum_error(solV3.q, solP3.q[end]) == 0

# end


# @testset "$(rpad("VPRK integrators with projection on Runge-Kutta tableau",80))" begin

#     # TODO: reactivate

#     # int = IntegratorVPRKpTableau(iode, CoefficientsPGLRK(5), Δt*5)
#     # sol = integrate(iode, int, div(nt,5))
#     # @test relative_maximum_error(sol.q, refx) < 8E-12

#     # int = IntegratorVPRKpTableau(iode, CoefficientsPGLRK(6), Δt*5)
#     # sol = integrate(iode, int, div(nt,5))
#     # @test relative_maximum_error(sol.q, refx) < 4E-14

# end


# @testset "$(rpad("Degenerate symplectic partitioned Runge-Kutta methods",80))" begin

#     sol = integrate(iode, DegenerateVPRK(VPRKGauss(1)))
#     @test relative_maximum_error(sol.q, reference_solution) < 2E-5

#     sol = integrate(iode, DegenerateVPRK(VPRKGauss(2)))
#     @test relative_maximum_error(sol.q, reference_solution) < 4E-7

#     sol = integrate(iode, DegenerateVPRK(VPRKGauss(3)))
#     @test relative_maximum_error(sol.q, reference_solution) < 2E-10

# end


# @testset "$(rpad("VPRK integrators with Legendre projection",80))" begin

#     sol = integrate(iode, VPRKpLegendre(VPRKGauss(1)))
#     @test relative_maximum_error(sol.q, reference_solution) < 1E-6

#     sol = integrate(iode, VPRKpLegendre(VPRKGauss(2)))
#     @test relative_maximum_error(sol.q, reference_solution) < 1E-11

#     sol = integrate(iode, VPRKpLegendre(VPRKGauss(3)))
#     @test relative_maximum_error(sol.q, reference_solution) < 4E-15

# end

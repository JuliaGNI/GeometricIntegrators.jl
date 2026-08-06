using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using Test

const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt * nt)

ode = odeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
iode = iodeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)
lode = lodeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)

ref = integrate(ode, Gauss(8))


@testset "$(rpad("Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, Gauss(1))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(iode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(iode, Gauss(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-11

    sol = integrate(iode, Gauss(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-13


    sol = integrate(iode, IRK(Gauss(1); implicit_update=true))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(iode, IRK(Gauss(2); implicit_update=true))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(iode, IRK(Gauss(3); implicit_update=true))
    @test relative_maximum_error(sol.q, ref.q) < 4E-11

    sol = integrate(iode, IRK(Gauss(4); implicit_update=true))
    @test relative_maximum_error(sol.q, ref.q) < 2E-13

end


@testset "$(rpad("Partitioned Runge-Kutta integrators for implicit equations",80))" begin

    sol = integrate(iode, PartitionedGauss(1))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(iode, PartitionedGauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(iode, PartitionedGauss(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-11

    sol = integrate(iode, PartitionedGauss(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-13

end


@testset "$(rpad("Formal Lagrangian Runge-Kutta integrators",80))" begin

    # The position phase of FLRK is exactly an implicit Runge-Kutta method applied to
    # q̇ = v̄(q), so it must reproduce `integrate(ode, Gauss(s))` to round-off. This is a
    # much stronger check than the error bound below, and it holds for every s.
    # (Not bit-for-bit in general: FLRK seeds Newton from `v̄` on top of the Hermite
    # extrapolation, so the converged iterates may differ by an ulp at coarse Δt. At
    # this Δt the results happen to be identical.)
    for s in 2:4
        @test relative_maximum_error(integrate(lode, FLRK(Gauss(s))).q,
                                     integrate(ode, Gauss(s)).q) < 4E-16
    end

    # Tolerances measured (Δt = 0.01, nt = 10, reference Gauss(8)):
    #   s = 2: 3.3E-12,  s = 3: 6.7E-16,  s = 4: 5.5E-16
    flsol = integrate(lode, FLRK(Gauss(2)))
    @test relative_maximum_error(flsol.q, ref.q) < 4E-12

    flsol = integrate(lode, FLRK(Gauss(3)))
    @test relative_maximum_error(flsol.q, ref.q) < 8E-16

    flsol = integrate(lode, FLRK(Gauss(4)))
    @test relative_maximum_error(flsol.q, ref.q) < 8E-16

    # The momentum-modified formal Lagrangian is constructed so that p = ϑ(q) solves
    # the adjoint equations. With p₀ = ϑ(q₀) the adjoint variable must therefore track
    # ϑ(q) along the whole trajectory. Measured: 1.0E-11 / 4.4E-16 / 1.6E-15.
    #
    # This assertion is what catches the stage force being taken from the projection
    # field `g = (∇ϑ)ᵀv` instead of `f = (∇ϑ)ᵀv − ∇H`: with `g` the error is 9.9E-2 at
    # every order, while `q` stays bit-for-bit correct.
    for (s, tol) in ((2, 2E-11), (3, 8E-16), (4, 2E-15))
        flsol = integrate(lode, FLRK(Gauss(s)))
        ϑerr = maximum(maximum(abs, flsol.p[i] .- LotkaVolterra2d.ϑ(flsol.t[i], flsol.q[i]))
                       for i in eachindex(flsol.q))
        @test ϑerr < tol
    end

end

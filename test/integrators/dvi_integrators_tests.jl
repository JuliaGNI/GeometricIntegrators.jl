using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dSingular as LVSingular
using SimpleSolvers
using Test
using Logging
using LinearAlgebra: I


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


@testset "$(rpad("1st Order DVIs",80))" begin

    sol = integrate(lode, DVIA())
    @test relative_maximum_error(sol.q, ref.q) < 1E-1

    sol = integrate(lode, DVIB())
    @test relative_maximum_error(sol.q, ref.q) < 1E-1

end


@testset "$(rpad("2nd Order Centred DVIs",80))" begin

    sol = integrate(lode, CMDVI())
    @test relative_maximum_error(sol.q, ref.q) < 4E-3

    sol = integrate(lode, CTDVI())
    @test relative_maximum_error(sol.q, ref.q) < 4E-3

end


@testset "$(rpad("Degenerate Variational Runge-Kutta integrators",80))" begin

    # NOTE: `LotkaVolterra2d` uses the gauge ϑ = (q₂ + log q₂/q₁, q₁), in which
    # neither component of ϑ vanishes. This is *outside* the class DVRK is built
    # for; the method is convergent but only of order s. See the in-class testset
    # below and `docs/src/integrators/dvi.md`.
    sol = integrate(lode, DVRK(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-5

    sol = integrate(lode, DVRK(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-7

    sol = integrate(lode, DVRK(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-10

    sol = integrate(lode, DVRK(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-13

end


@testset "$(rpad("DVRK on an in-class degenerate Lagrangian",80))" begin

    # `LotkaVolterra2dSingular` uses the gauge ϑ = (log q₂/q₁, 0), which satisfies
    # ϑ_μ = 0 for μ > d/2. This is the class for which DVRK is symplectic and of
    # full order 2s. Same dynamics as `LotkaVolterra2d` above, hence same reference.
    slode = LVSingular.lodeproblem(q₀; timespan=tspan, timestep=Δt, parameters=params)

    sol = integrate(slode, DVRK(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(slode, DVRK(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(slode, DVRK(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-14

end


@testset "$(rpad("DVRK preserves the noncanonical symplectic form",80))" begin

    # The defining property of DVRK: the one-step map q_n ↦ q_{n+1} preserves
    # ω = dϑ, i.e. Jᵀ Ω(q₁) J = Ω(q₀) with J the Jacobian of the map and
    # Ω_μν = ∂ϑ_ν/∂q^μ - ∂ϑ_μ/∂q^ν. This holds only under the hypotheses of the
    # theorem, two of which are exercised below: ϑ in a gauge with ϑ_μ = 0 for
    # μ > d/2, and a tableau satisfying b_i a_ij + b_j a_ji = b_i b_j.
    Ω(m, q) = (o = m.dϑ₂dx₁(0.0, q) - m.dϑ₁dx₂(0.0, q); [0.0 -o; o 0.0])

    function onestep(m, method, q, h)
        prob = m.lodeproblem(collect(q); timespan=(0.0, h), timestep=h, parameters=params)
        collect(integrate(prob, method).q[end])
    end

    function symplecticity_defect(m, method, q₀, h; δ=1e-6)
        f(q) = onestep(m, method, q, h)
        J = reduce(hcat, [(f(q₀ .+ δ .* e) .- f(q₀ .- δ .* e)) ./ 2δ
                          for e in eachcol(Matrix(1.0I, 2, 2))])
        maximum(abs, J' * Ω(m, f(q₀)) * J .- Ω(m, q₀))
    end

    qref = [2.0, 1.5]

    # In class with a symplectic tableau: exact to the finite-difference floor,
    # and — crucially — independent of the step size.
    for h in (0.1, 0.5), s in 1:3
        @test symplecticity_defect(LVSingular, DVRK(Gauss(s)), qref, h) < 1e-8
    end

    # Out of class (ϑ₂ = q₁ ≠ 0): the same dynamics in a different gauge, no
    # longer symplectic, and the defect grows with the step size. Note that the
    # invertibility hypothesis is *not* what fails here — ∂ϑ₁/∂q² = 1 + 1/q₁q₂ ≠ 0
    # for this gauge too; what fails is ϑ_μ = 0 for μ > d/2.
    @test symplecticity_defect(LotkaVolterra2d, DVRK(Gauss(2)), qref, 0.1) > 1e-6
    @test symplecticity_defect(LotkaVolterra2d, DVRK(Gauss(2)), qref, 0.5) >
          symplecticity_defect(LotkaVolterra2d, DVRK(Gauss(2)), qref, 0.1)

    # In class, but with a tableau violating b_i a_ij + b_j a_ji = b_i b_j.
    @test symplecticity_defect(LVSingular, DVRK(RadauIIA(2); check_conditions=false), qref, 0.5) > 1e-6

end


@testset "$(rpad("DVRK method traits and tableau validation",80))" begin

    # Regression: these used to call an undefined `tableaus(method)` and threw
    # `UndefVarError` on any DVRK instance.
    @test order(DVRK(Gauss(1))) == 2
    @test order(DVRK(Gauss(2))) == 4
    @test order(DVRK(Gauss(3))) == 6
    @test issymmetric(DVRK(Gauss(2)))
    @test issymplectic(DVRK(Gauss(2)))
    @test isimplicit(DVRK(Gauss(2)))
    @test !isexplicit(DVRK(Gauss(2)))

    # Symplecticity of DVRK depends on the tableau, so unlike the other DVI
    # methods it cannot be answered for the bare type.
    @test ismissing(issymplectic(DVRK))
    @test issymplectic(CMDVI)

    # Tableaus violating the theorem's hypotheses are accepted but warn.
    @test !issymplectic(DVRK(RadauIIA(3); check_conditions=false))
    @test_logs (:warn, r"symplecticity condition") DVRK(RadauIIA(3))
    @test_logs (:warn, r"symplecticity condition") (:warn, r"singular") DVRK(LobattoIIIA(3))
    @test_nowarn DVRK(Gauss(2))

    # Regression: the invertibility check must not be `det`-based. det(a) of the
    # Gauss tableaus falls below eps^(3/4) from s = 10 on, although cond(a) ≈ 1e2.
    @test_nowarn DVRK(Gauss(10))
    @test_nowarn DVRK(Gauss(12))

    # Inconsistent initial momentum p₀ ≠ ϑ(q₀) warns; the default p₀ does not.
    @test_logs min_level=Logging.Warn LVSingular.lodeproblem([2.0, 3.0]; timespan=tspan, timestep=Δt) |>
        prob -> integrate(prob, DVRK(Gauss(1)))
    @test_logs (:warn, r"not consistent with the primary constraint") match_mode=:any begin
        integrate(LVSingular.lodeproblem([2.0, 3.0], [0.0, 0.0]; timespan=tspan, timestep=Δt), DVRK(Gauss(1)))
    end

end

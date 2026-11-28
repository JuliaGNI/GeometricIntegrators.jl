using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

include("../helpers/verification_utilities.jl")

# The harmonic oscillator has an analytic solution, so we can measure the true
# error directly. It is a *linear* problem, so it measures the linear order, and
# some methods super-converge on it. The harmonic-oscillator testsets below
# therefore only require the empirical order to match the documented order
# (super-convergence beyond that is not flagged); only the testset "Linear order
# of Lobatto IIIF, IIIF̄ and IIIG" asserts a linear order above the classical one.
const T = 1.0
build(Δt) = odeproblem(; timespan = (0.0, T), timestep = Δt)

# Refinement levels: coarser base steps for higher-order methods so that the
# coarsest error stays above the roundoff plateau (the plateau filter in
# `estimate_convergence_order` discards the saturated tail).
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))

# A nonlinear, non-autonomous, non-stiff scalar problem with a closed-form solution:
# y' = cos(t)⋅y², y(0) = 1/2, exact y(t) = 1/(2 - sin t) ∈ [1/2, 1] and |∂f/∂y| ≤ 2
# on [0, T].
# On it a method shows its classical order, which the harmonic oscillator cannot
# separate from the linear order, and a wrong stage time shows as a lost order.
nonlinear_v(v, t, q, params) = (v[1] = cos(t) * q[1]^2; nothing)
nonlinear_build(Δt) = ODEProblem(nonlinear_v, (0.0, T), Δt, [0.5])
nonlinear_exact(t) = 1 / (2 - sin(t))

# Maximum over all output times of |q - y(t)|, relative to the maximum of |y(t)|.
function nonlinear_error(sol, ref)
    err = maximum(i -> abs(sol.q[i][1] - nonlinear_exact(sol.t[i])), eachindex(sol.t))
    return err / maximum(i -> abs(nonlinear_exact(sol.t[i])), eachindex(sol.t))
end

@testset "Runge-Kutta convergence" begin
    @testset "Explicit Runge-Kutta" begin
        test_convergence_order(
            build, ExplicitEulerRK(), steps(10, 5); reference = exact_solution,
            expected = 1, label = "ExplicitEulerRK")
        test_convergence_order(
            build, ExplicitMidpoint(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "ExplicitMidpoint")
        test_convergence_order(build, Heun2(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "Heun2")
        test_convergence_order(build, Heun3(), steps(10, 5); reference = exact_solution,
            expected = 3, label = "Heun3")
        test_convergence_order(build, Kutta3(), steps(10, 5); reference = exact_solution,
            expected = 3, label = "Kutta3")
        test_convergence_order(build, Ralston2(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "Ralston2")
        test_convergence_order(build, Ralston3(), steps(10, 5); reference = exact_solution,
            expected = 3, label = "Ralston3")
        test_convergence_order(build, Runge2(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "Runge2")
        test_convergence_order(build, SSPRK2(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "SSPRK2")
        test_convergence_order(build, SSPRK3(), steps(10, 5); reference = exact_solution,
            expected = 3, label = "SSPRK3")
        test_convergence_order(build, RK4(), steps(10, 4); reference = exact_solution,
            expected = 4, label = "RK4")
        test_convergence_order(build, RK438(), steps(10, 4); reference = exact_solution,
            expected = 4, label = "RK438")
    end

    @testset "Diagonally implicit Runge-Kutta" begin
        test_convergence_order(
            build, CrankNicolsonRK(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "CrankNicolsonRK")
        test_convergence_order(build, Crouzeix(), steps(10, 5); reference = exact_solution,
            expected = 3, label = "Crouzeix")
        test_convergence_order(build, QinZhang(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "QinZhang")
        # KraaijevangerSpijker is a first-order method (its order attribute was
        # corrected to 1 in RungeKutta.jl v0.5.22).
        test_convergence_order(
            build, KraaijevangerSpijker(), steps(10, 5); reference = exact_solution,
            expected = 1, label = "KraaijevangerSpijker")
    end

    @testset "Fully implicit Runge-Kutta" begin
        test_convergence_order(
            build, ImplicitEulerRK(), steps(10, 5); reference = exact_solution,
            expected = 1, label = "ImplicitEulerRK")
        test_convergence_order(
            build, ImplicitMidpointRK(), steps(10, 5); reference = exact_solution,
            expected = 2, label = "ImplicitMidpointRK")
        test_convergence_order(build, SRK3(), steps(4, 4); reference = exact_solution,
            expected = 4, label = "SRK3")
        test_convergence_order(build, Gauss(1), steps(10, 4); reference = exact_solution,
            expected = 2, label = "Gauss(1)")
        test_convergence_order(build, Gauss(2), steps(4, 4); reference = exact_solution,
            expected = 4, label = "Gauss(2)")
        test_convergence_order(build, Gauss(3), steps(2, 4); reference = exact_solution,
            expected = 6, label = "Gauss(3)")
    end

    @testset "Radau and Lobatto" begin
        test_convergence_order(build, RadauIA(2), steps(6, 4); reference = exact_solution,
            expected = 3, label = "RadauIA(2)")
        test_convergence_order(build, RadauIIA(2), steps(6, 4); reference = exact_solution,
            expected = 3, label = "RadauIIA(2)")
        test_convergence_order(
            build, LobattoIIIA(2), steps(10, 4); reference = exact_solution,
            expected = 2, label = "LobattoIIIA(2)")
        test_convergence_order(
            build, LobattoIIIC(2), steps(10, 4); reference = exact_solution,
            expected = 2, label = "LobattoIIIC(2)")
        test_convergence_order(
            build, LobattoIIID(2), steps(10, 4); reference = exact_solution,
            expected = 2, label = "LobattoIIID(2)")
        test_convergence_order(
            build, LobattoIIIE(2), steps(10, 4); reference = exact_solution,
            expected = 2, label = "LobattoIIIE(2)")
        test_convergence_order(
            build, LobattoIIIA(3), steps(6, 4); reference = exact_solution,
            expected = 4, label = "LobattoIIIA(3)")
        test_convergence_order(
            build, LobattoIIIB(3), steps(6, 4); reference = exact_solution,
            expected = 4, label = "LobattoIIIB(3)")
    end

    # Lobatto IIIF and IIIF̄ have the (s, s)-Padé stability function, so on the linear
    # harmonic oscillator they reach the linear order 2s, above their classical order
    # 2s - 2; Lobatto IIIG has linear order 2s - 2.
    @testset "Linear order of Lobatto IIIF, IIIF̄ and IIIG" begin
        test_convergence_order(
            build, LobattoIIIF(2), steps(10, 4); reference = exact_solution,
            expected = 4, label = "LobattoIIIF(2), linear order")
        test_convergence_order(
            build, LobattoIIIF̄(2), steps(10, 4); reference = exact_solution,
            expected = 4, label = "LobattoIIIF̄(2), linear order")
        test_convergence_order(
            build, LobattoIIIG(2), steps(10, 4); reference = exact_solution,
            expected = 2, label = "LobattoIIIG(2), linear order")
    end

    # Classical order on the nonlinear non-autonomous problem, with the same base
    # steps and refinement levels as on the harmonic oscillator.
    @testset "Nonlinear non-autonomous problem" begin
        cases = (
            ("ExplicitEulerRK", ExplicitEulerRK(), 10, 5, 1),
            ("ExplicitMidpoint", ExplicitMidpoint(), 10, 5, 2),
            ("Heun2", Heun2(), 10, 5, 2),
            ("Heun3", Heun3(), 10, 5, 3),
            ("Kutta3", Kutta3(), 10, 5, 3),
            ("Ralston2", Ralston2(), 10, 5, 2),
            ("Ralston3", Ralston3(), 10, 5, 3),
            ("Runge2", Runge2(), 10, 5, 2),
            ("SSPRK2", SSPRK2(), 10, 5, 2),
            ("SSPRK3", SSPRK3(), 10, 5, 3),
            ("RK4", RK4(), 10, 4, 4),
            ("RK438", RK438(), 10, 4, 4),
            ("CrankNicolsonRK", CrankNicolsonRK(), 10, 5, 2),
            ("Crouzeix", Crouzeix(), 10, 5, 3),
            ("QinZhang", QinZhang(), 10, 5, 2),
            ("KraaijevangerSpijker", KraaijevangerSpijker(), 10, 5, 1),
            ("ImplicitEulerRK", ImplicitEulerRK(), 10, 5, 1),
            ("ImplicitMidpointRK", ImplicitMidpointRK(), 10, 5, 2),
            ("SRK3", SRK3(), 4, 4, 4),
            ("Gauss(1)", Gauss(1), 10, 4, 2),
            ("Gauss(2)", Gauss(2), 4, 4, 4),
            ("Gauss(3)", Gauss(3), 2, 4, 6),
            ("RadauIA(2)", RadauIA(2), 6, 4, 3),
            ("RadauIIA(2)", RadauIIA(2), 6, 4, 3),
            ("LobattoIIIA(2)", LobattoIIIA(2), 10, 4, 2),
            ("LobattoIIIC(2)", LobattoIIIC(2), 10, 4, 2),
            ("LobattoIIID(2)", LobattoIIID(2), 10, 4, 2),
            ("LobattoIIIE(2)", LobattoIIIE(2), 10, 4, 2),
            ("LobattoIIIA(3)", LobattoIIIA(3), 6, 4, 4),
            ("LobattoIIIB(3)", LobattoIIIB(3), 6, 4, 4),
            ("LobattoIIIF(2)", LobattoIIIF(2), 10, 4, 2),
            ("LobattoIIIF(3)", LobattoIIIF(3), 6, 4, 4),
            ("LobattoIIIF̄(2)", LobattoIIIF̄(2), 10, 4, 2),
            ("LobattoIIIF̄(3)", LobattoIIIF̄(3), 6, 4, 4),
            ("LobattoIIIG(2)", LobattoIIIG(2), 10, 4, 2),
            ("LobattoIIIG(3)", LobattoIIIG(3), 6, 4, 4)
        )
        for (label, method, n0, k, expected) in cases
            test_convergence_order(
                nonlinear_build, method, steps(n0, k); reference = prob -> nothing,
                errormetric = nonlinear_error, expected = expected,
                label = label * ", nonlinear")
        end
    end

    # The classical order of the Lobatto IIIF, IIIF̄ and IIIG tableaus and of their
    # partitioned and VPRK forms is 2s - 2: their Lobatto weights satisfy B(2s - 2)
    # only. The instances take their order from the RungeKutta tableau.
    @testset "Order of the Lobatto IIIF, IIIF̄ and IIIG methods" begin
        @testset "order($(nameof(M)))" for M in (
            LobattoIIIF, LobattoIIIF̄, LobattoIIIG,
            LobattoIIIFIIIF̄, LobattoIIIF̄IIIF, LobattoIIIGIIIḠ,
            VPRKLobattoIIIF, VPRKLobattoIIIF̄, VPRKLobattoIIIG,
            VPRKLobattoIIIFIIIF̄, VPRKLobattoIIIF̄IIIF, VPRKLobattoIIIGIIIḠ)
            @test order(M) == "2s-2"
        end
        @testset "order($(nameof(M))($s))" for M in (
                LobattoIIIF, LobattoIIIF̄, LobattoIIIG,
                LobattoIIIFIIIF̄, LobattoIIIF̄IIIF, LobattoIIIGIIIḠ),
            s in 2:4

            @test order(M(s)) == 2s - 2
        end
    end

    # Known order deficiency (audit finding 3 in docs/src/audit.md): LobattoIIIB
    # has a singular coefficient matrix, so as a standalone ODE integrator its stage
    # system is degenerate and it drops to order 1, on the harmonic oscillator and on
    # the nonlinear problem alike (it is intended for use in symplectic partitioned
    # pairs / as a VPRK component). Recorded as broken so a future change is
    # detected automatically.
    @testset "Known order deficiencies (broken)" begin
        r = estimate_convergence_order(build, LobattoIIIB(2), steps(10, 4); reference = exact_solution)
        @test_broken isapprox(r.order, 2; atol = 0.35)  # issue #247
    end
end

using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

include("verification_utilities.jl")

# The harmonic oscillator has an analytic solution, so we can measure the true
# error directly. It is a *linear* problem, so some methods super-converge on it;
# the assertions below therefore only require the empirical order to match the
# documented order (super-convergence beyond that is not flagged).
const T = 1.0
build(Δt) = odeproblem(; timespan = (0.0, T), timestep = Δt)

# Refinement levels: coarser base steps for higher-order methods so that the
# coarsest error stays above the roundoff plateau (the plateau filter in
# `estimate_convergence_order` discards the saturated tail).
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))

@testset "Runge-Kutta convergence" begin

    @testset "Explicit Runge-Kutta" begin
        test_convergence_order(build, ExplicitEulerRK(), steps(10, 5); reference = exact_solution, expected = 1, label = "ExplicitEulerRK")
        test_convergence_order(build, ExplicitMidpoint(), steps(10, 5); reference = exact_solution, expected = 2, label = "ExplicitMidpoint")
        test_convergence_order(build, Heun2(),  steps(10, 5); reference = exact_solution, expected = 2, label = "Heun2")
        test_convergence_order(build, Heun3(),  steps(10, 5); reference = exact_solution, expected = 3, label = "Heun3")
        test_convergence_order(build, Kutta3(), steps(10, 5); reference = exact_solution, expected = 3, label = "Kutta3")
        test_convergence_order(build, Ralston2(), steps(10, 5); reference = exact_solution, expected = 2, label = "Ralston2")
        test_convergence_order(build, Ralston3(), steps(10, 5); reference = exact_solution, expected = 3, label = "Ralston3")
        test_convergence_order(build, Runge2(), steps(10, 5); reference = exact_solution, expected = 2, label = "Runge2")
        test_convergence_order(build, SSPRK2(), steps(10, 5); reference = exact_solution, expected = 2, label = "SSPRK2")
        test_convergence_order(build, SSPRK3(), steps(10, 5); reference = exact_solution, expected = 3, label = "SSPRK3")
        test_convergence_order(build, RK4(),   steps(10, 4); reference = exact_solution, expected = 4, label = "RK4")
        test_convergence_order(build, RK438(), steps(10, 4); reference = exact_solution, expected = 4, label = "RK438")
    end

    @testset "Diagonally implicit Runge-Kutta" begin
        test_convergence_order(build, CrankNicolson(), steps(10, 5); reference = exact_solution, expected = 2, label = "CrankNicolson")
        test_convergence_order(build, Crouzeix(),      steps(10, 5); reference = exact_solution, expected = 3, label = "Crouzeix")
        test_convergence_order(build, QinZhang(),      steps(10, 5); reference = exact_solution, expected = 2, label = "QinZhang")
    end

    @testset "Fully implicit Runge-Kutta" begin
        test_convergence_order(build, ImplicitEulerRK(), steps(10, 5); reference = exact_solution, expected = 1, label = "ImplicitEulerRK")
        test_convergence_order(build, ImplicitMidpoint(), steps(10, 5); reference = exact_solution, expected = 2, label = "ImplicitMidpoint")
        test_convergence_order(build, SRK3(), steps(4, 4); reference = exact_solution, expected = 4, label = "SRK3")
        test_convergence_order(build, Gauss(1), steps(10, 4); reference = exact_solution, expected = 2, label = "Gauss(1)")
        test_convergence_order(build, Gauss(2), steps(4, 4);  reference = exact_solution, expected = 4, label = "Gauss(2)")
        test_convergence_order(build, Gauss(3), steps(2, 4);  reference = exact_solution, expected = 6, label = "Gauss(3)")
    end

    @testset "Radau and Lobatto" begin
        test_convergence_order(build, RadauIA(2),  steps(6, 4); reference = exact_solution, expected = 3, label = "RadauIA(2)")
        test_convergence_order(build, RadauIIA(2), steps(6, 4); reference = exact_solution, expected = 3, label = "RadauIIA(2)")
        test_convergence_order(build, LobattoIIIA(2), steps(10, 4); reference = exact_solution, expected = 2, label = "LobattoIIIA(2)")
        test_convergence_order(build, LobattoIIIC(2), steps(10, 4); reference = exact_solution, expected = 2, label = "LobattoIIIC(2)")
        test_convergence_order(build, LobattoIIID(2), steps(10, 4); reference = exact_solution, expected = 2, label = "LobattoIIID(2)")
        test_convergence_order(build, LobattoIIIE(2), steps(10, 4); reference = exact_solution, expected = 2, label = "LobattoIIIE(2)")
        test_convergence_order(build, LobattoIIIA(3), steps(6, 4); reference = exact_solution, expected = 4, label = "LobattoIIIA(3)")
        test_convergence_order(build, LobattoIIIB(3), steps(6, 4); reference = exact_solution, expected = 4, label = "LobattoIIIB(3)")
    end

    # Known order deficiencies (see VERIFICATION_REPORT.md). These methods do not
    # reach their documented order when applied as standalone ODE integrators;
    # the root cause is in the tableau coefficients (RungeKutta.jl). Recorded as
    # broken so a future fix is detected automatically.
    @testset "Known order deficiencies (broken)" begin
        # KraaijevangerSpijker: coefficients satisfy only the order-1 conditions
        # (Σᵢbᵢcᵢ = 2 ≠ 1/2), yet order(·) and the docs claim order 2.
        r = estimate_convergence_order(build, KraaijevangerSpijker(), steps(10, 5); reference = exact_solution)
        @test_broken isapprox(r.order, 2; atol = 0.35)
        # LobattoIIIB(2): singular A ⇒ degenerate stage system ⇒ order 1 standalone.
        r = estimate_convergence_order(build, LobattoIIIB(2), steps(10, 4); reference = exact_solution)
        @test_broken isapprox(r.order, 2; atol = 0.35)
    end

end

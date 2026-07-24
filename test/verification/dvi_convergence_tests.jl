using GeometricIntegrators
import GeometricProblems.Pendulum as Pendulum
using GeometricProblems.LotkaVolterra2d
using Test

include("verification_utilities.jl")

const T = 1.0
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

@testset "Degenerate variational integrator convergence" begin

    # DVRK on the (regular) pendulum IODE reaches the full order 2s of its
    # underlying Gauss-Legendre tableau.
    @testset "DVRK order (regular Lagrangian)" begin
        vbuild(Δt) = Pendulum.iodeproblem(; timespan = (0.0, T), timestep = Δt)
        vref(prob) = integrate(prob, VPRKGauss(8))
        test_convergence_order(vbuild, DVRK(Gauss(2)), steps(4, 3); reference = vref, errormetric = emq, expected = 4, label = "DVRK(Gauss(2))")
        test_convergence_order(vbuild, DVRK(Gauss(3)), steps(2, 3); reference = vref, errormetric = emq, expected = 6, label = "DVRK(Gauss(3))")
    end

    # On the degenerate Lotka-Volterra Lagrangian, DVRK exhibits order reduction:
    # DVRK(Gauss(s)) converges at order s rather than 2s (a known feature of
    # symplectic RK methods on degenerate Lagrangians; see VERIFICATION_REPORT.md).
    @testset "DVRK order reduction (degenerate Lagrangian)" begin
        q₀ = [1.0, 1.0]
        params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
        lvbuild(Δt) = lodeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
        lvref(prob) = integrate(odeproblem(q₀; timespan = timespan(prob), timestep = timestep(prob), parameters = params), Gauss(8))
        test_convergence_order(lvbuild, DVRK(Gauss(1)), steps(20, 4); reference = lvref, errormetric = emq, expected = 1, label = "DVRK(Gauss(1))/LV")
        test_convergence_order(lvbuild, DVRK(Gauss(2)), steps(10, 4); reference = lvref, errormetric = emq, expected = 2, plateau = 1e-11, label = "DVRK(Gauss(2))/LV")
    end

    # Fixed-step accuracy regression for the low-order degenerate variational
    # integrators (matching the established Lotka-Volterra test). NOTE: over long
    # integration times these methods do not converge cleanly to the ODE reference;
    # see VERIFICATION_REPORT.md. Here we only assert the short-time accuracy that
    # the package already relies on.
    @testset "Low-order DVI accuracy (fixed step)" begin
        q₀ = [1.0, 1.0]
        params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
        tspan = (0.0, 0.1)
        lode = lodeproblem(q₀; timespan = tspan, timestep = 0.01, parameters = params)
        ref = integrate(odeproblem(q₀; timespan = tspan, timestep = 0.01, parameters = params), Gauss(8))
        @test relative_maximum_error(integrate(lode, DVIA()).q, ref.q) < 1e-1
        @test relative_maximum_error(integrate(lode, DVIB()).q, ref.q) < 1e-1
        @test relative_maximum_error(integrate(lode, CMDVI()).q, ref.q) < 4e-3
        @test relative_maximum_error(integrate(lode, CTDVI()).q, ref.q) < 4e-3
    end
end

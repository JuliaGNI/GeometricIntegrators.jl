using GeometricIntegrators
import GeometricProblems.Pendulum as Pendulum
import GeometricProblems.LotkaVolterra2dSingular as LVSingular
import GeometricProblems.LotkaVolterra2dSymmetric as LVSymmetric
import GeometricProblems.MasslessChargedParticle as MCP
import GeometricProblems.MasslessChargedParticleSingular as MCPSingular
using GeometricProblems.LotkaVolterra2d
using Test

include("verification_utilities.jl")

const T = 1.0
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

@testset "Degenerate variational integrator convergence" begin

    # The pendulum IODE is the degenerate phase-space form of the pendulum: on
    # q = (x, p) ∈ R² it has ϑ = (m l² q₂, 0), so ϑ₂ ≡ 0 and ∂ϑ₁/∂q₂ = m l² ≠ 0.
    # It is therefore *in* the DVRK class, and reaches the full order 2s of its
    # underlying Gauss-Legendre tableau.
    @testset "DVRK order (degenerate pendulum, in class)" begin
        vbuild(Δt) = Pendulum.iodeproblem(; timespan = (0.0, T), timestep = Δt)
        vref(prob) = integrate(prob, VPRKGauss(8))
        test_convergence_order(vbuild, DVRK(Gauss(2)), steps(4, 3); reference = vref, errormetric = emq, expected = 4, label = "DVRK(Gauss(2))")
        test_convergence_order(vbuild, DVRK(Gauss(3)), steps(2, 3); reference = vref, errormetric = emq, expected = 6, label = "DVRK(Gauss(3))")
    end

    # DVRK is designed for degenerate Lagrangians L = ϑ(q)⋅q̇ - H(q) in which d/2
    # components of ϑ vanish identically. For such Lagrangians it attains the full
    # order 2s of the underlying Gauss-Legendre tableau. `LotkaVolterra2dSingular`
    # (ϑ₂ = 0) and `MasslessChargedParticleSingular` (A₂ = 0) are in this class.
    @testset "DVRK order (degenerate Lagrangian, in class)" begin
        q₀ = [2.0, 1.0]
        params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
        lvbuild(Δt) = LVSingular.lodeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
        lvref(prob) = integrate(odeproblem(q₀; timespan = timespan(prob), timestep = timestep(prob), parameters = params), Gauss(8))
        test_convergence_order(lvbuild, DVRK(Gauss(1)), steps(20, 4); reference = lvref, errormetric = emq, expected = 2, label = "DVRK(Gauss(1))/LV-singular")
        test_convergence_order(lvbuild, DVRK(Gauss(2)), steps(10, 4); reference = lvref, errormetric = emq, expected = 4, label = "DVRK(Gauss(2))/LV-singular")
        test_convergence_order(lvbuild, DVRK(Gauss(3)), steps(5, 4); reference = lvref, errormetric = emq, expected = 6, plateau = 1e-12, label = "DVRK(Gauss(3))/LV-singular")

        mcpbuild(Δt) = MCPSingular.lodeproblem(; timespan = (0.0, T), timestep = Δt)
        mcpref(prob) = integrate(MCP.odeproblem(; timespan = timespan(prob), timestep = timestep(prob)), Gauss(8))
        test_convergence_order(mcpbuild, DVRK(Gauss(1)), steps(20, 4); reference = mcpref, errormetric = emq, expected = 2, label = "DVRK(Gauss(1))/MCP-singular")
        test_convergence_order(mcpbuild, DVRK(Gauss(2)), steps(10, 4); reference = mcpref, errormetric = emq, expected = 4, plateau = 1e-12, label = "DVRK(Gauss(2))/MCP-singular")
    end

    # The symplectic potential ϑ is only defined up to an exact one-form, and the
    # gauge matters: written in a gauge in which no component of ϑ vanishes, the
    # *same* dynamics falls outside the class DVRK is built for and the order drops
    # from 2s to s. These problems describe exactly the same systems as above.
    @testset "DVRK order reduction (degenerate Lagrangian, out of class)" begin
        q₀ = [2.0, 1.0]
        params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
        lvref(prob) = integrate(odeproblem(q₀; timespan = timespan(prob), timestep = timestep(prob), parameters = params), Gauss(8))

        # ϑ = (q₂ + log q₂/q₁, q₁): gauge-equivalent to the singular form, ϑ₂ ≠ 0
        lvbuild(Δt) = lodeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
        test_convergence_order(lvbuild, DVRK(Gauss(1)), steps(20, 4); reference = lvref, errormetric = emq, expected = 1, label = "DVRK(Gauss(1))/LV")
        test_convergence_order(lvbuild, DVRK(Gauss(2)), steps(10, 4); reference = lvref, errormetric = emq, expected = 2, label = "DVRK(Gauss(2))/LV")

        # ϑ = (log q₂/2q₁, -log q₁/2q₂): also gauge-equivalent, both components ≠ 0
        symbuild(Δt) = LVSymmetric.lodeproblem(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
        test_convergence_order(symbuild, DVRK(Gauss(2)), steps(10, 4); reference = lvref, errormetric = emq, expected = 2, label = "DVRK(Gauss(2))/LV-symmetric")

        # A₁, A₂ both ≠ 0: gauge-equivalent to MasslessChargedParticleSingular
        mcpbuild(Δt) = MCP.lodeproblem(; timespan = (0.0, T), timestep = Δt)
        mcpref(prob) = integrate(MCP.odeproblem(; timespan = timespan(prob), timestep = timestep(prob)), Gauss(8))
        test_convergence_order(mcpbuild, DVRK(Gauss(2)), steps(10, 4); reference = mcpref, errormetric = emq, expected = 2, label = "DVRK(Gauss(2))/MCP")
    end

    # Fixed-step accuracy regression for the low-order degenerate variational
    # integrators (matching the established Lotka-Volterra test). NOTE: over long
    # integration times these methods do not converge cleanly to the ODE reference;
    # see docs/src/audit.md. Here we only assert the short-time accuracy that
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

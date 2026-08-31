using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
import GeometricProblems.Pendulum as Pendulum
using Test

include("verification_utilities.jl")

const T = 1.0
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

# Position-momentum and discrete-Euler-Lagrange integrators on the harmonic
# oscillator LODE, referenced against the analytic PODE solution.
lbuild(Δt) = lodeproblem(; timespan = (0.0, T), timestep = Δt)
function pref(prob)
    exact_solution(podeproblem(; timespan = timespan(prob), timestep = timestep(prob)))
end
delbuild_m(Δt) = DELEProblem(lodeproblem(; timespan = (0.0, T), timestep = Δt), Midpoint())
function delbuild_t(Δt)
    DELEProblem(lodeproblem(; timespan = (0.0, T), timestep = Δt), Trapezoidal())
end

# VPRK integrators on the nonlinear pendulum IODE (true nonlinear order). The high-order
# VPRKGauss(8) reference is relaxed to f_abstol = 4e-15 for the same reason as PMVI below:
# at the finest timestep its residual stagnates at the round-off floor, which otherwise
# trips one spurious stagnation warning.
vbuild(Δt) = Pendulum.iodeproblem(; timespan = (0.0, T), timestep = Δt)
vref(prob) = integrate(prob, VPRKGauss(8); f_abstol = 4e-15)

@testset "Variational integrator convergence" begin
    @testset "Position-momentum and discrete Euler-Lagrange" begin
        # PMVI relaxes the solver residual tolerance to f_abstol = 4e-15 to avoid a
        # spurious stagnation warning at the finest timestep, where the residual settles
        # at a round-off floor of ~4e-15 — above the framework default
        # f_abstol = max(8, solversize)·eps, here its floor 8eps ≈ 1.78e-15 (see
        # docs/src/audit.md). The DEL solves stay below that floor and need no relaxation.
        # @test_nowarn is the tripwire: a converged solve must not warn.
        @test_nowarn test_convergence_order(
            lbuild, PMVImidpoint(), steps(10, 4); reference = pref,
            errormetric = emq, expected = 2, label = "PMVImidpoint",
            integrate_options = (f_abstol = 4e-15,))
        @test_nowarn test_convergence_order(
            lbuild, PMVItrapezoidal(), steps(10, 4); reference = pref,
            errormetric = emq, expected = 2, label = "PMVItrapezoidal",
            integrate_options = (f_abstol = 4e-15,))
        @test_nowarn test_convergence_order(
            delbuild_m, DiscreteEulerLagrange(), steps(10, 4); reference = pref,
            errormetric = emq, expected = 2, label = "DEL-midpoint")
        @test_nowarn test_convergence_order(
            delbuild_t, DiscreteEulerLagrange(), steps(10, 4); reference = pref,
            errormetric = emq, expected = 2, label = "DEL-trapezoidal")
    end

    @testset "VPRK" begin
        # VPRKGauss(2) reaches the round-off floor at the finest of its steps(5, 4)
        # timesteps; f_abstol = 4e-15 (as for PMVI above) keeps the converged solve silent
        # instead of tripping a stagnation warning.
        @test_nowarn test_convergence_order(
            vbuild, VPRKGauss(2), steps(5, 4); reference = vref,
            errormetric = emq, expected = 4, label = "VPRKGauss(2)",
            integrate_options = (f_abstol = 4e-15,))
        @test_nowarn test_convergence_order(
            vbuild, VPRKGauss(3), steps(3, 3); reference = vref,
            errormetric = emq, expected = 6, label = "VPRKGauss(3)")
        @test_nowarn test_convergence_order(
            vbuild, VPRKLobattoIIIAIIIB(2), steps(5, 4); reference = vref,
            errormetric = emq, expected = 2, label = "VPRKLobattoIIIAIIIB(2)")

        # Known order deficiencies (see docs/src/audit.md): several higher-order
        # Lobatto VPRK methods do not reach their documented order. The "order 2s"
        # IIIF/IIIG methods only reach order 2 for s = 2, and VPRKLobattoIIIAIIIB(3)
        # reaches order 2 rather than 2s-2 = 4 (whereas the plain partitioned
        # LobattoIIIAIIIB(3) does reach order 4).
        @testset "Known order deficiencies (broken)" begin
            r = estimate_convergence_order(vbuild, VPRKLobattoIIIF̄(2), steps(5, 4);
                reference = vref, errormetric = emq)
            @test_broken isapprox(r.order, 4; atol = 0.4)
            r = estimate_convergence_order(vbuild, VPRKLobattoIIIG(2), steps(5, 4);
                reference = vref, errormetric = emq)
            @test_broken isapprox(r.order, 4; atol = 0.4)
            r = estimate_convergence_order(vbuild, VPRKLobattoIIIAIIIB(3), steps(3, 4);
                reference = vref, errormetric = emq)
            @test_broken isapprox(r.order, 4; atol = 0.4)
        end
    end
end

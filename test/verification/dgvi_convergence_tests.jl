using GeometricIntegrators
import GeometricProblems.LotkaVolterra2d as LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dGauge as LotkaVolterra2dGauge
using CompactBasisFunctions
using QuadratureRules
using Test

include("verification_utilities.jl")

# Discontinuous Galerkin variational integrators on the degenerate Lotka-Volterra
# Lagrangian, referenced against a high-order solution of the equivalent ODE.
#
# The five variants split cleanly into two groups by how they discretise the numerical
# flux, and the measured orders differ sharply:
#
#   * `DGVIPI` (path integral) and `DGVIEXP` (midpoint average) evaluate the one-form at an
#     *average* of the two one-sided limits and reach the full order 2s.
#   * `DGVI`, `DGVIP0` and `DGVIP1` evaluate a trapezoidal flux involving the *nodal*
#     value qₙ and are capped at 2⌊s/2⌋ — order 2 for s = 2 and s = 3, order 4 for s = 4.
#
# Measured (T = 1, Δt = 1/4 … 1/32):
#
#     s | DGVI  DGVIP0  DGVIP1 | DGVIPI  DGVIEXP
#     2 | 1.94   2.01    1.94  |  1.97    2.01
#     3 | 2.01   2.01    2.01  |  5.95    5.97
#     4 | 3.74   3.86    3.74  |  5.97*   6.00*      (* machine-precision limited)
#
# The deficiency of the nodal-value flux is consistent with the manuscript's own
# `% TODO The jump condition is not quite right` annotation. It is recorded here as
# `@test_broken` against the nominal 2s rather than asserted away.

const T = 1.0
const q₀ = [1.0, 1.0]
const params = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)

function build(Δt)
    LotkaVolterra2d.iodeproblem_dg(q₀; timespan = (0.0, T), timestep = Δt, parameters = params)
end
gbuild(Δt) = LotkaVolterra2dGauge.iodeproblem_dg(; timespan = (0.0, T), timestep = Δt)
function ref(prob)
    integrate(
        LotkaVolterra2d.odeproblem(
            q₀; timespan = timespan(prob), timestep = timestep(prob), parameters = params),
        Gauss(8))
end
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, r) = relative_maximum_error(sol.q, r.q)

basisquad(s) = (Q = GaussLegendreQuadrature(s); (Lagrange(QuadratureRules.nodes(Q)), Q))
const JUMP = Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))

dgvi(s) = DGVI(basisquad(s)...)
dgvipi(s) = DGVIPI(basisquad(s)..., JUMP)
dgvip0(s) = DGVIP0(basisquad(s)...)
dgvip1(s) = DGVIP1(basisquad(s)...)
dgviexp(s) = DGVIEXP(basisquad(s)...)

@testset "$(rpad("Discontinuous Galerkin variational integrator convergence",80))" begin
    @testset "average-based flux reaches order 2s" begin
        test_convergence_order(build, dgvipi(2), steps(4, 3); reference = ref,
            errormetric = emq, expected = 2, label = "DGVIPI(2)")
        test_convergence_order(build, dgvipi(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 6, label = "DGVIPI(3)")

        test_convergence_order(build, dgviexp(2), steps(4, 3); reference = ref,
            errormetric = emq, expected = 2, label = "DGVIEXP(2)")
        test_convergence_order(build, dgviexp(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 6, label = "DGVIEXP(3)")
    end

    @testset "nodal-value flux is capped at 2⌊s/2⌋" begin
        # what these variants actually achieve
        test_convergence_order(build, dgvi(2), steps(4, 3); reference = ref,
            errormetric = emq, expected = 2, label = "DGVI(2)")
        test_convergence_order(build, dgvi(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 2, label = "DGVI(3)")
        test_convergence_order(build, dgvi(4), steps(4, 3); reference = ref,
            errormetric = emq, expected = 4, atol = 0.4, label = "DGVI(4)")

        test_convergence_order(build, dgvip0(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 2, label = "DGVIP0(3)")
        test_convergence_order(build, dgvip1(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 2, label = "DGVIP1(3)")

        # what the nominal order would be: the flux involving the nodal value qₙ is only
        # second-order accurate and limits the whole scheme.
        for (label, m) in (("DGVI(3)", dgvi(3)), ("DGVIP0(3)", dgvip0(3)), (
            "DGVIP1(3)", dgvip1(3)))
            r = estimate_convergence_order(
                build, m, steps(4, 3); reference = ref, errormetric = emq)
            @test_broken isapprox(r.order, 6; atol = 0.35)
        end
    end

    @testset "gauge-transformed problem" begin
        # The Lagrangian may be augmented by the total time derivative of a gauge term
        # without changing the continuous Euler-Lagrange equations, and a variational
        # integrator must be insensitive to that. `iodeproblem_dg_gauge` carries the gauge
        # parameter κ for exactly this purpose; κ = 0 recovers `iodeproblem_dg`.
        κbuild(κ) = Δt -> LotkaVolterra2d.iodeproblem_dg_gauge(
            q₀; timespan = (0.0, T), timestep = Δt, parameters = params, κ = κ)

        test_convergence_order(κbuild(0.0), dgviexp(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 6, atol = 0.5, label = "DGVIEXP(3) κ=0")
        test_convergence_order(κbuild(0.5), dgviexp(3), steps(4, 3); reference = ref,
            errormetric = emq, expected = 6, atol = 0.5, label = "DGVIEXP(3) κ=1/2")

        # `LotkaVolterra2dGauge` writes the same system with a gauge-transformed one-form
        # in a separate module; it must agree too.
        gref(prob) = integrate(
            LotkaVolterra2dGauge.odeproblem(; timespan = timespan(prob), timestep = timestep(prob)),
            Gauss(8))
        test_convergence_order(
            gbuild, dgviexp(3), steps(4, 3); reference = gref, errormetric = emq,
            expected = 6, atol = 0.5, label = "DGVIEXP(3) gauge module")
    end
end

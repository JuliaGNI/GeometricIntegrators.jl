using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using QuadratureRules
using CompactBasisFunctions
using Test

include("verification_utilities.jl")

# Continuous Galerkin variational integrators with a Lagrange basis, on the harmonic
# oscillator IODE (analytic PODE reference).
#
# `CGVI` is run on Gauss-Legendre nodes and `CGVINodal` on Lobatto-Legendre nodes — the
# latter requires the interval endpoints, which Gauss nodes do not include. Both converge
# at order 2s-2 in the position, so the two families are asserted against the same numbers.
const T = 1.0
build(Δt) = iodeproblem(; timespan = (0.0, T), timestep = Δt)
pref(prob) = exact_solution(podeproblem(; timespan = timespan(prob), timestep = timestep(prob)))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

cgvi(s) = (Q = GaussLegendreQuadrature(s); CGVI(Lagrange(QuadratureRules.nodes(Q)), Q))
cgvinodal(s) = (Q = LobattoLegendreQuadrature(s); CGVINodal(Lagrange(QuadratureRules.nodes(Q)), Q))

@testset "Continuous Galerkin variational integrator convergence" begin
    test_convergence_order(build, cgvi(2), steps(2, 4); reference = pref, errormetric = emq, expected = 2, label = "CGVI(Gauss(2))")
    test_convergence_order(build, cgvi(3), steps(2, 3); reference = pref, errormetric = emq, expected = 4, label = "CGVI(Gauss(3))")
    test_convergence_order(build, cgvi(4), steps(2, 3); reference = pref, errormetric = emq, expected = 6, label = "CGVI(Gauss(4))")

    test_convergence_order(build, cgvinodal(2), steps(2, 4); reference = pref, errormetric = emq, expected = 2, label = "CGVINodal(Lobatto(2))")
    test_convergence_order(build, cgvinodal(3), steps(2, 3); reference = pref, errormetric = emq, expected = 4, label = "CGVINodal(Lobatto(3))")
    test_convergence_order(build, cgvinodal(4), steps(2, 3); reference = pref, errormetric = emq, expected = 6, label = "CGVINodal(Lobatto(4))")
end

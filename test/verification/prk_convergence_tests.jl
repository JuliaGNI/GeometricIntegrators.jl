using GeometricIntegrators
import GeometricProblems.Pendulum as Pendulum
using Test

include("verification_utilities.jl")

# Nonlinear pendulum PODE with a high-order reference, so that the measured order
# is the true (nonlinear) order rather than a linear super-convergence artifact.
const T = 1.0
build(Δt) = Pendulum.podeproblem(; timespan = (0.0, T), timestep = Δt)
href(prob) = integrate(prob, PartitionedGauss(8))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

@testset "Partitioned Runge-Kutta convergence" begin
    test_convergence_order(build, SymplecticEulerA(), steps(20, 4); reference = href, errormetric = emq, expected = 1, label = "SymplecticEulerA")
    test_convergence_order(build, SymplecticEulerB(), steps(20, 4); reference = href, errormetric = emq, expected = 1, label = "SymplecticEulerB")
    test_convergence_order(build, PartitionedGauss(2), steps(5, 4); reference = href, errormetric = emq, expected = 4, label = "PartitionedGauss(2)")
    test_convergence_order(build, LobattoIIIAIIIB(2), steps(5, 4); reference = href, errormetric = emq, expected = 2, label = "LobattoIIIAIIIB(2)")
    test_convergence_order(build, LobattoIIIBIIIA(2), steps(5, 4); reference = href, errormetric = emq, expected = 2, label = "LobattoIIIBIIIA(2)")
    test_convergence_order(build, LobattoIIIAIIIB(3), steps(5, 3); reference = href, errormetric = emq, expected = 4, label = "LobattoIIIAIIIB(3)")

    # Known order deficiencies (see VERIFICATION_REPORT.md): the "order 2s" Lobatto
    # IIIF/IIIG pairs only reach order 2s-2 = 2 for s = 2.
    @testset "Known order deficiencies (broken)" begin
        r = estimate_convergence_order(build, LobattoIIIFIIIF̄(2), steps(5, 4); reference = href, errormetric = emq)
        @test_broken isapprox(r.order, 4; atol = 0.4)
        r = estimate_convergence_order(build, LobattoIIIGIIIḠ(2), steps(5, 4); reference = href, errormetric = emq)
        @test_broken isapprox(r.order, 4; atol = 0.4)
    end
end

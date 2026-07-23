using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

include("verification_utilities.jl")

# Splitting/composition methods are tested on the harmonic oscillator SODE; the
# reference is the analytic solution of the corresponding ODE.
const T = 1.0
build(Δt) = sodeproblem(; timespan = (0.0, T), timestep = Δt)
sref(prob) = exact_solution(odeproblem(; timespan = timespan(prob), timestep = timestep(prob)))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

@testset "Splitting and composition convergence" begin
    test_convergence_order(build, LieA(),   steps(10, 4); reference = sref, errormetric = emq, expected = 1, label = "LieA")
    test_convergence_order(build, LieB(),   steps(10, 4); reference = sref, errormetric = emq, expected = 1, label = "LieB")
    test_convergence_order(build, Strang(), steps(10, 4); reference = sref, errormetric = emq, expected = 2, label = "Strang")
    test_convergence_order(build, McLachlan2(), steps(10, 4); reference = sref, errormetric = emq, expected = 2, label = "McLachlan2")
    test_convergence_order(build, McLachlan4(), steps(10, 4); reference = sref, errormetric = emq, expected = 4, label = "McLachlan4")
    test_convergence_order(build, TripleJump(), steps(10, 4); reference = sref, errormetric = emq, expected = 4, label = "TripleJump")
    test_convergence_order(build, SuzukiFractal(), steps(10, 4); reference = sref, errormetric = emq, expected = 4, label = "SuzukiFractal")

    # Higher-order Yoshida composition methods (reach the roundoff plateau quickly,
    # so use coarser step ranges / a longer interval).
    test_convergence_order(build, Yoshida6(), steps(4, 4); reference = sref, errormetric = emq, expected = 6, label = "Yoshida6")
    let T2 = 2.0
        build2(Δt) = sodeproblem(; timespan = (0.0, T2), timestep = Δt)
        sref2(prob) = exact_solution(odeproblem(; timespan = timespan(prob), timestep = timestep(prob)))
        test_convergence_order(build2, Yoshida8(), T2 ./ (2 .* 2 .^ (0:5)); reference = sref2, errormetric = emq, expected = 8, label = "Yoshida8")
    end
end

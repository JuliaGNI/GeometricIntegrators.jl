using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

include("verification_utilities.jl")

const T = 1.0
build(Δt) = lodeproblem(; timespan = (0.0, T), timestep = Δt)
pref(prob) = exact_solution(podeproblem(; timespan = timespan(prob), timestep = timestep(prob)))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, ref) = relative_maximum_error(sol.q, ref.q)

# Trivial velocity map φ_h(q̄, q) = (q - q̄)/Δt and its derivatives, matching the
# existing Hamilton-Pontryagin integrator tests.
ϕ(v, q̄, q, a, Δt) = (v .= (q .- q̄) ./ Δt)
function D₁ϕ(d, q̄, q, a, Δt)
    d .= 0
    for i in eachindex(q)
        d[i, i] = -1 / Δt
    end
end
function D₂ϕ(d, q̄, q, a, Δt)
    d .= 0
    for i in eachindex(q)
        d[i, i] = +1 / Δt
    end
end
Dₐϕ(d, q̄, q, a, Δt) = (d .= 0)

@testset "Hamilton-Pontryagin convergence" begin
    # The HP integrators relax the solver residual tolerance to f_abstol = 4e-15 to
    # avoid a spurious "Solver took 1000 iterations." warning at the finest timestep
    # (tolerance stagnation near machine precision; see VERIFICATION_REPORT.md, third
    # pass). min_iterations = 1 is repeated because any solver option replaces the
    # whole default_options bundle.
    test_convergence_order(build, HPImidpoint(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[]), steps(10, 4);
        reference = pref, errormetric = emq, expected = 2, label = "HPImidpoint",
        integrate_options = (min_iterations = 1, f_abstol = 4e-15))
    test_convergence_order(build, HPItrapezoidal(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[]), steps(10, 4);
        reference = pref, errormetric = emq, expected = 2, label = "HPItrapezoidal",
        integrate_options = (min_iterations = 1, f_abstol = 4e-15))
end

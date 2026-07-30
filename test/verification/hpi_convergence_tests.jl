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
    # avoid a spurious stagnation warning at the finest timestep, where the residual
    # settles at the round-off floor above the default (effectively zero) f_abstol
    # (see docs/src/audit.md, third pass). GeometricIntegratorsBase merges this into
    # the method's default_options, so only f_abstol needs to be given. @test_nowarn
    # is a regression tripwire: under SimpleSolvers 0.10 these converged solves are
    # silent, so any resurfacing solver warning fails here.
    @test_nowarn test_convergence_order(build, HPImidpoint(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[]), steps(10, 4);
        reference = pref, errormetric = emq, expected = 2, label = "HPImidpoint",
        integrate_options = (f_abstol = 4e-15,))
    @test_nowarn test_convergence_order(build, HPItrapezoidal(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[]), steps(10, 4);
        reference = pref, errormetric = emq, expected = 2, label = "HPItrapezoidal",
        integrate_options = (f_abstol = 4e-15,))
end

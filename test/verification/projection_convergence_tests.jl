using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

include("verification_utilities.jl")

# Projected Runge-Kutta integrators on the harmonic oscillator DAE, referenced
# against the analytic ODE solution. Projected Gauss(s) converges at order 2s and
# preserves the Hamiltonian to machine precision.
const T = 1.0
build(Δt) = daeproblem(; timespan = (0.0, T), timestep = Δt)
ref(prob) = exact_solution(odeproblem(; timespan = timespan(prob), timestep = timestep(prob)))
steps(n0, k) = T ./ (n0 .* 2 .^ (0:k))
emq(sol, r) = relative_maximum_error(sol.q, r.q)

energy_change(sol, params) = abs(hamiltonian(sol[end].t, sol[end].q, params) -
                                 hamiltonian(sol[begin].t, sol[begin].q, params))

@testset "Projected integrator convergence and energy conservation" begin
    for (P, name, etol) in ((PostProjection, "PostProjection", 1),
                            (MidpointProjection, "MidpointProjection", 2),
                            (SymmetricProjection, "SymmetricProjection", 2))
        @testset "$name" begin
            test_convergence_order(build, P(Gauss(1)), steps(10, 4); reference = ref, errormetric = emq, expected = 2, label = "$name(Gauss(1))")
            test_convergence_order(build, P(Gauss(2)), steps(5, 4);  reference = ref, errormetric = emq, expected = 4, label = "$name(Gauss(2))")
            test_convergence_order(build, P(Gauss(3)), steps(3, 3);  reference = ref, errormetric = emq, expected = 6, label = "$name(Gauss(3))")

            # Energy conservation over the trajectory (to machine precision).
            prob = build(0.1)
            for s in 1:3
                sol = integrate(prob, P(Gauss(s)))
                @test energy_change(sol, parameters(prob)) < etol * eps()
            end
        end
    end
end

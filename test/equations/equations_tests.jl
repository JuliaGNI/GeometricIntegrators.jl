
using GeometricIntegrators.Equations
using Test

t₀ = 0.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]


include("deterministic_equations_tests.jl")
include("stochastic_equations_tests.jl")

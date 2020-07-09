
using GeometricIntegrators.Common
using GeometricIntegrators.Equations
using GeometricIntegrators.Equations: function_v_dummy, get_λ₀
using Test

include("initial_conditions.jl")

zero_vec = [zeros(3) for i in 1:3]

@test function_v_dummy(t₀, q₀, p₀, λ₀) == nothing

@test get_λ₀(zeros(3), zeros(3)) == zeros(3)
@test get_λ₀(zeros(3), zero_vec) == zeros(3)
@test get_λ₀(zero_vec, zeros(3)) == zero_vec
@test get_λ₀(zero_vec, zero_vec) == zero_vec

struct TestEquation{DT,TT} <: Equation{DT,TT} end

@test_throws ErrorException ndims(TestEquation{Float64,Float64}())

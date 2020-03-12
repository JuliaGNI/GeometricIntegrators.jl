
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Equations
using GeometricIntegrators.Equations: function_v_dummy, get_λ₀
using Test

t₀ = 0.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]

qₛ = rand(1,3,3)
pₛ = rand(1,3,3)
xₛ = rand(2,3,3)


@testset "$(rpad("General equation functionality",80))" begin

    @test function_v_dummy(t₀, q₀, p₀, λ₀) == nothing

    @test get_λ₀(zeros(3), zeros(3))   == zeros(3)
    @test get_λ₀(zeros(3), zeros(3,3)) == zeros(3)
    @test get_λ₀(zeros(3,3), zeros(3)) == zeros(3,3)

    struct TestEquation{DT,TT} <: Equation{DT,TT} end

    @test_throws ErrorException ndims(TestEquation{Float64,Float64}())

end


include("deterministic_equations_tests.jl")
include("stochastic_equations_tests.jl")

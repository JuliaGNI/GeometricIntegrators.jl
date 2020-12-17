
using GeometricIntegrators.Common
using GeometricIntegrators.Equations
using Test

include("initial_conditions.jl")


@testset "$(rpad("Stochastic Differential Equations (SDE)",80))" begin

    x₁ₛ = [0.5  0.0 -0.5; 0.0  0.5  0.0]

    function sde_v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function sde_B(μ, t, q, u)
        u[1] = μ*q[1]
        u[2] = μ*q[2]
    end

    λ  = 2.
    μ  = 1.

    sde_v_params = (t, q, v) -> sde_v(λ, t, q, v)
    sde_B_params = (t, q, v) -> sde_B(μ, t, q, v)

    sde  = SDE(1, 1, sde_v_params, sde_B_params, t₀, x₀)
    sde1 = SDE(1, 1, sde_v_params, sde_B_params, x₀)
    sde2 = SDE(1, 3, sde_v_params, sde_B_params, x₀)
    sde3 = SDE(1, 1, sde_v_params, sde_B_params, x₁ₛ)

    @test ndims(sde) == 2
    @test periodicity(sde) == zero(x₀)
    @test get_function_tuple(sde) == NamedTuple{(:v,:B)}((sde_v_params, sde_B_params))

    @test sde == sde1
    @test sde != sde2
    @test sde != sde3

    @test hash(sde) == hash(sde1)
    @test hash(sde) != hash(sde2)
    @test hash(sde) != hash(sde3)

    @test sde1.d == 2
    @test sde2.d == 2
    @test sde3.d == 2

    @test sde1.ni == 1
    @test sde2.ni == 1
    @test sde3.ni == 3

    @test sde1.ns == 1
    @test sde2.ns == 3
    @test sde3.ns == 1

    @test sde == similar(sde, t₀, x₀, 1)
    @test sde == similar(sde, t₀, x₀)
    @test sde == similar(sde, x₀, 1)
    @test sde == similar(sde, x₀)

end


@testset "$(rpad("Partitioned Stochastic Differential Equations (PSDE)",80))" begin

    noise_intensity = 0.1

    q₁ₛ = [0.5  0.0 -0.5]
    p₁ₛ = [0.0  0.5  0.0]

    function psde_v(t, q, p, v)
        v[1] = p[1]
    end

    function psde_f(t, q, p, f)
        f[1] = - q[1]
    end

    function psde_B(t, q, p, B)
        B[1,1] = noise_intensity * p[1]
    end

    function psde_G(t, q, p, G)
        G[1,1] = - noise_intensity * q[1]
    end

    psde  = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, t₀, q₀, p₀)
    psde1 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, q₀, p₀)
    psde2 = PSDE(1, 3, psde_v, psde_f, psde_B, psde_G, q₀, p₀)
    psde3 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, q₁ₛ, p₁ₛ)

    @test ndims(psde) == 1
    @test periodicity(psde) == zero(q₀)
    @test get_function_tuple(psde) == NamedTuple{(:v,:f,:B,:G)}((psde_v,psde_f,psde_B,psde_G))

    @test psde == psde1
    @test psde != psde2
    @test psde != psde3

    @test hash(psde) == hash(psde1)
    @test hash(psde) != hash(psde2)
    @test hash(psde) != hash(psde3)

    @test psde1.d == 1
    @test psde2.d == 1
    @test psde3.d == 1

    @test psde1.ni == 1
    @test psde2.ni == 1
    @test psde3.ni == 3

    @test psde1.ns == 1
    @test psde2.ns == 3
    @test psde3.ns == 1

    @test psde == similar(psde, t₀, q₀, p₀, 1)
    @test psde == similar(psde, q₀, p₀, 1)

    @test psde != similar(psde, t₀, q₁ₛ, p₁ₛ)
    @test psde != similar(psde, q₁ₛ, p₁ₛ)

end


@testset "$(rpad("Split Partitioned Stochastic Differential Equations (SPSDE)",80))" begin

    noise_intensity = 0.1

    q₁ₛ = [0.5  0.0 -0.5]
    p₁ₛ = [0.0  0.5  0.0]

    function spsde_v(t, q, p, v)
        v[1] = p[1]
    end

    function spsde_f1(t, q, p, f)
        f[1] = - q[1]
    end

    function spsde_f2(t, q, p, f)
        f[1] = 0
    end

    function spsde_B(t, q, p, B)
        B[1,1] = noise_intensity * p[1]
    end

    function spsde_G1(t, q, p, G)
        G[1,1] = - noise_intensity * q[1]
    end

    function spsde_G2(t, q, p, G)
        G[1,1] = 0
    end

    spsde  = SPSDE(1, 1, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, t₀, q₀, p₀)
    spsde1 = SPSDE(1, 1, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, q₀, p₀)
    spsde2 = SPSDE(1, 3, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, q₀, p₀)
    spsde3 = SPSDE(1, 1, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, q₁ₛ, p₁ₛ)

    @test ndims(spsde) == 1
    @test periodicity(spsde) == zero(q₀)
    @test get_function_tuple(spsde) == NamedTuple{(:v,:f1,:f2,:B,:G1,:G2)}((spsde_v,spsde_f1,spsde_f2,spsde_B,spsde_G1,spsde_G2))

    @test spsde == spsde1
    @test spsde != spsde2
    @test spsde != spsde3

    @test hash(spsde) == hash(spsde1)
    @test hash(spsde) != hash(spsde2)
    @test hash(spsde) != hash(spsde3)

    @test spsde1.d == 1
    @test spsde2.d == 1
    @test spsde3.d == 1

    @test spsde1.ni == 1
    @test spsde2.ni == 1
    @test spsde3.ni == 3

    @test spsde1.ns == 1
    @test spsde2.ns == 3
    @test spsde3.ns == 1

    @test spsde == similar(spsde, t₀, q₀, p₀, 1)
    @test spsde == similar(spsde, q₀, p₀, 1)

    @test spsde != similar(spsde, t₀, q₁ₛ, p₁ₛ)
    @test spsde != similar(spsde, q₁ₛ, p₁ₛ)

end

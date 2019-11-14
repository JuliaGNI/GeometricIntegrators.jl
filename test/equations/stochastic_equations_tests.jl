
@testset "$(rpad("Stochastic equations",80))" begin

    ################################################################################
    # Test SDE: Stochastic Differential Equation
    ################################################################################

    function v_sde(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function u_sde(μ, t, q, u)
        u[1] = μ*q[1]
        u[2] = μ*q[2]
    end

    λ  = 2.
    μ  = 1.

    v_sde_params = (t, q, v) -> v_sde(λ, t, q, v)
    u_sde_params = (t, q, v) -> u_sde(μ, t, q, v)

    sde  = SDE{eltype(x₀), typeof(t₀), typeof(v_sde_params), typeof(u_sde_params), 1}(2, 1, 1, 1, v_sde_params, u_sde_params, t₀, x₀)
    sde1 = SDE(1, 1, v_sde_params, u_sde_params, t₀, x₀)
    sde2 = SDE(1, 1, v_sde_params, u_sde_params, x₀)

    @test sde == sde1
    @test sde == sde2

    @test hash(sde1) == hash(sde2)

    @test sde == similar(sde, t₀, x₀, 1)
    @test sde == similar(sde, t₀, x₀)
    @test sde == similar(sde, x₀, 1)
    @test sde == similar(sde, x₀)

end

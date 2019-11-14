
@testset "$(rpad("Stochastic equations",80))" begin

    ################################################################################
    # Test SDE: Stochastic Differential Equation
    ################################################################################

    function sde_v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function sde_u(μ, t, q, u)
        u[1] = μ*q[1]
        u[2] = μ*q[2]
    end

    λ  = 2.
    μ  = 1.

    sde_v_params = (t, q, v) -> sde_v(λ, t, q, v)
    sde_u_params = (t, q, v) -> sde_u(μ, t, q, v)

    sde  = SDE{eltype(x₀), typeof(t₀), typeof(sde_v_params), typeof(sde_u_params), 1}(2, 1, 1, 1, sde_v_params, sde_u_params, t₀, x₀)
    sde1 = SDE(1, 1, sde_v_params, sde_u_params, t₀, x₀)
    sde2 = SDE(1, 1, sde_v_params, sde_u_params, x₀)

    @test sde == sde1
    @test sde == sde2

    @test hash(sde1) == hash(sde2)

    @test sde == similar(sde, t₀, x₀, 1)
    @test sde == similar(sde, t₀, x₀)
    @test sde == similar(sde, x₀, 1)
    @test sde == similar(sde, x₀)

end

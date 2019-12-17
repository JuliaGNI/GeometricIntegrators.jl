
@testset "$(rpad("Stochastic equations",80))" begin

    ################################################################################
    # Test SDE: Stochastic Differential Equation
    ################################################################################

    x₁ₛ = [0.5  0.0 -0.5; 0.0  0.5  0.0]

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

    sde  = SDE(1, 1, 1, sde_v_params, sde_u_params, t₀, x₀)
    sde1 = SDE(1, 1, sde_v_params, sde_u_params, t₀, x₀)
    sde2 = SDE(1, 1, sde_v_params, sde_u_params, x₀)
    sde3 = SDE(1, 1, sde_v_params, sde_u_params, x₁ₛ)
    sde4 = SDE(1, sde_v_params, sde_u_params, x₁ₛ)
    sde5 = SDE(1, sde_v_params, sde_u_params, xₛ)

    @test sde == sde1
    @test sde == sde2
    @test sde != sde3
    @test sde != sde4
    @test sde != sde5

    @test hash(sde1) == hash(sde2)
    @test hash(sde3) != hash(sde4)

    @test sde1.d == 2
    @test sde2.d == 2
    @test sde3.d == 2
    @test sde4.d == 2
    @test sde5.d == 2

    @test sde1.n == 1
    @test sde2.n == 1
    @test sde3.n == 3
    @test sde4.n == 1
    @test sde5.n == 3

    @test sde1.ns == 1
    @test sde2.ns == 1
    @test sde3.ns == 1
    @test sde4.ns == 3
    @test sde5.ns == 3

    @test sde == similar(sde, t₀, x₀, 1)
    @test sde == similar(sde, t₀, x₀)
    @test sde == similar(sde, x₀, 1)
    @test sde == similar(sde, x₀)


    ################################################################################
    # Test PSDE: Partitioned Stochastic Differential Equation
    ################################################################################

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

    psde  = PSDE(1, 1, 1, psde_v, psde_f, psde_B, psde_G, t₀, q₀, p₀)
    psde1 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, t₀, q₀, p₀)
    psde2 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, q₀, p₀)
    psde3 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, q₁ₛ, p₁ₛ)
    psde4 = PSDE(1, psde_v, psde_f, psde_B, psde_G, q₁ₛ, p₁ₛ)
    psde5 = PSDE(1, psde_v, psde_f, psde_B, psde_G, qₛ, pₛ)

    @test psde == psde1
    @test psde == psde2
    @test psde != psde3
    @test psde != psde4
    @test psde != psde5

    @test hash(psde1) == hash(psde2)
    @test hash(psde3) != hash(psde4)

    @test psde1.d == 1
    @test psde2.d == 1
    @test psde3.d == 1
    @test psde4.d == 1
    @test psde5.d == 1

    @test psde1.n == 1
    @test psde2.n == 1
    @test psde3.n == 3
    @test psde4.n == 1
    @test psde5.n == 3

    @test psde1.ns == 1
    @test psde2.ns == 1
    @test psde3.ns == 1
    @test psde4.ns == 3
    @test psde5.ns == 3

    @test psde == similar(psde, t₀, q₀, p₀, 1)
    @test psde == similar(psde, q₀, p₀, 1)

    @test psde != similar(psde, t₀, q₁ₛ, p₁ₛ)
    @test psde != similar(psde, q₁ₛ, p₁ₛ)

end

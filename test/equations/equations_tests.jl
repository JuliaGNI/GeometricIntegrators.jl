
t₀ = 0.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]

@testset "$(rpad("Deterministic equations",80))" begin

    ################################################################################
    # Test ODE: Ordinary Differential Equation
    ################################################################################

    function f_ode(t, x, f)
        f[1] = x[1]
    end

    ode  = ODE(eltype(q₀), 1, 1, 1, f_ode, t₀, q₀)
    ode1 = ODE(f_ode, t₀, q₀)
    ode2 = ODE(f_ode, q₀)

    @test ode == ode1
    @test ode == ode2

    @test hash(ode1) == hash(ode2)

    @test ode == similar(ode, t₀, q₀)
    @test ode == similar(ode, q₀)


    ################################################################################
    # Test PODE: Partitioned Ordinary Differential Equation
    ################################################################################

    function v_pode(t, q, p, v)
        v[1] = q[1]
    end

    function f_pode(t, q, p, f)
        f[1] = 2p[1]
    end

    pode  = PODE(eltype(q₀), 1, 1, 1, v_pode, f_pode, t₀, q₀, p₀)
    pode1 = PODE(v_pode, f_pode, t₀, q₀, p₀)
    pode2 = PODE(v_pode, f_pode, q₀, p₀)

    @test pode == pode1
    @test pode == pode2

    @test hash(pode1) == hash(pode2)

    @test pode == similar(pode, t₀, q₀, p₀)
    @test pode == similar(pode, q₀, p₀)


    ################################################################################
    # Test IODE: Implicit Ordinary Differential Equation
    ################################################################################

    function iode_ϑ(t, q, v, p)
        p[1] = v[1]
    end

    function iode_f(t, q, v, f)
        f[1] = sin(q[1])
    end

    function iode_g(t, q, λ, g)
        g[1] = λ[1]
    end

    function iode_v(t, q, p, v)
        v[1] = p[1]
    end

    iode  = IODE(eltype(q₀), 1, 1, 1, iode_ϑ, iode_f, iode_g, t₀, q₀, p₀, λ₀; v=iode_v)
    iode1 = IODE(iode_ϑ, iode_f, iode_g, t₀, q₀, p₀, λ₀; v=iode_v)
    iode2 = IODE(iode_ϑ, iode_f, iode_g, t₀, q₀, p₀; v=iode_v)
    iode3 = IODE(iode_ϑ, iode_f, iode_g, q₀, p₀; v=iode_v)

    @test iode == iode1
    @test iode == iode2
    @test iode == iode3

    @test hash(iode1) == hash(iode2)

    @test iode == similar(iode, t₀, q₀, p₀, λ₀)
    @test iode == similar(iode, t₀, q₀, p₀)
    @test iode == similar(iode, q₀, p₀)


    ################################################################################
    # Test DAE: Differential Algebraic Equation
    ################################################################################

    function v_dae(t, x, v)
        v[1] = x[1]
        v[2] = x[2]
    end

    function u_dae(t, x, λ, u)
        u[1] = +λ[1]
        u[1] = -λ[1]
    end

    function ϕ_dae(t, x, λ, ϕ)
        ϕ[1] = x[2] - x[1]
    end

    dae  = DAE(eltype(q₀), 1, 2, 1, 1, v_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
    dae1 = DAE(v_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
    dae2 = DAE(v_dae, u_dae, ϕ_dae, x₀, λ₀)

    @test dae == dae1
    @test dae == dae2

    @test hash(dae1) == hash(dae2)

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)


    ################################################################################
    # Test PDAE: Partitioned Differential Algebraic Equation
    ################################################################################

    function v_pdae(t, q, p, v)
        v[1] = q[1]
    end

    function f_pdae(t, q, p, f)
        f[1] = p[1]
    end

    function p_pdae(t, q, v, p)
        p[1] = v[1]
    end

    function u_pdae(t, q, p, λ, u)
        u[1] = λ[1]
    end

    function g_pdae(t, q, p, λ, g)
        g[1] = λ[1]
    end

    function ϕ_pdae(t, q, p, λ, ϕ)
        ϕ[1] = p[1] - q[1]
    end

    pdae  = PDAE(eltype(q₀), 1, 1, 1, 1, v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
    pdae1 = PDAE(v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
    pdae2 = PDAE(v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, q₀, p₀, λ₀)

    @test pdae == pdae1
    @test pdae == pdae2

    @test hash(pdae1) == hash(pdae2)

    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)


    idae  = IDAE(eltype(q₀), 1, 1, 1, 1, p_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
    idae1 = IDAE(p_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
    idae2 = IDAE(p_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, q₀, p₀, λ₀)

    @test idae == idae1
    @test idae == idae2

    @test hash(idae1) == hash(idae2)

    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)

end

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

end

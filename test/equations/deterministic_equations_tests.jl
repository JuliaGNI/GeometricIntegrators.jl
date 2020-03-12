
@testset "$(rpad("Deterministic equations",80))" begin

    ################################################################################
    # Test ODE: Ordinary Differential Equation
    ################################################################################

    function ode_v(t, x, f)
        f[1] = x[1]
    end

    ode  = ODE(eltype(q₀), 1, 1, 1, ode_v, t₀, q₀)
    ode1 = ODE(ode_v, t₀, q₀)
    ode2 = ODE(ode_v, q₀)

    @test ode == ode1
    @test ode == ode2

    @test hash(ode1) == hash(ode2)

    @test ode == similar(ode, t₀, q₀)
    @test ode == similar(ode, q₀)


    ################################################################################
    # Test PODE: Partitioned Ordinary Differential Equation
    ################################################################################

    function pode_v(t, q, p, v)
        v[1] = q[1]
    end

    function pode_f(t, q, p, f)
        f[1] = 2p[1]
    end

    pode  = PODE(eltype(q₀), 1, 1, 1, pode_v, pode_f, t₀, q₀, p₀)
    pode1 = PODE(pode_v, pode_f, t₀, q₀, p₀)
    pode2 = PODE(pode_v, pode_f, q₀, p₀)

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
    # Test HODE: Hamiltonian Ordinary Differential Equation
    ################################################################################

    hode_eqs = (iode_v, iode_f, iode_h)

    hode  = HODE(eltype(q₀), 1, 1, 1, hode_eqs..., t₀, q₀, p₀)
    hode1 = HODE(hode_eqs..., t₀, q₀, p₀)
    hode2 = HODE(hode_eqs..., q₀, p₀)

    @test ndims(hode) == 1
    @test periodicity(hode) == zero(q₀)
    @test get_function_tuple(hode) == NamedTuple{(:v, :f, :h)}(hode_eqs)

    @test hode == hode1
    @test hode == hode2

    @test hash(hode1) == hash(hode2)

    @test hode == similar(hode, t₀, q₀, p₀)
    @test hode == similar(hode, q₀, p₀)


    ################################################################################
    # Test VODE: Variational Ordinary Differential Equation
    ################################################################################

    vode  = VODE(eltype(q₀), 1, 1, 1, iode_ϑ, iode_f, iode_g, t₀, q₀, p₀, λ₀; v=iode_v)
    vode1 = VODE(iode_ϑ, iode_f, iode_g, t₀, q₀, p₀, λ₀; v=iode_v)
    vode2 = VODE(iode_ϑ, iode_f, iode_g, t₀, q₀, p₀; v=iode_v)
    vode3 = VODE(iode_ϑ, iode_f, iode_g, q₀, p₀; v=iode_v)

    @test vode == vode1
    @test vode == vode2
    @test vode == vode3

    @test hash(vode1) == hash(vode2)

    @test vode == similar(vode, t₀, q₀, p₀, λ₀)
    @test vode == similar(vode, t₀, q₀, p₀)
    @test vode == similar(vode, q₀, p₀)


    ################################################################################
    # Test DAE: Differential Algebraic Equation
    ################################################################################

    function dae_v(t, x, v)
        v[1] = x[1]
        v[2] = x[2]
    end

    function dae_u(t, x, λ, u)
        u[1] = +λ[1]
        u[2] = -λ[1]
    end

    function dae_ϕ(t, x, λ, ϕ)
        ϕ[1] = x[2] - x[1]
    end

    dae  = DAE(eltype(q₀), 1, 2, 1, 1, dae_v, dae_u, dae_ϕ, t₀, x₀, λ₀)
    dae1 = DAE(dae_v, dae_u, dae_ϕ, t₀, x₀, λ₀)
    dae2 = DAE(dae_v, dae_u, dae_ϕ, x₀, λ₀)

    @test dae == dae1
    @test dae == dae2

    @test hash(dae1) == hash(dae2)

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)


    ################################################################################
    # Test PDAE: Partitioned Differential Algebraic Equation
    ################################################################################

    function pdae_v(t, q, p, v)
        v[1] = p[1]
    end

    function pdae_f(t, q, p, f)
        f[1] = q[1]
    end

    function pdae_p(t, q, v, p)
        p[1] = v[1]
    end

    function pdae_u(t, q, p, λ, u)
        u[1] = λ[1]
    end

    function pdae_g(t, q, p, λ, g)
        g[1] = λ[1]
    end

    function pdae_ϕ(t, q, p, λ, ϕ)
        ϕ[1] = p[1] - q[1]
    end

    pdae  = PDAE(eltype(q₀), 1, 1, 1, 1, pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, t₀, q₀, p₀, λ₀)
    pdae1 = PDAE(pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, t₀, q₀, p₀, λ₀)
    pdae2 = PDAE(pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, q₀, p₀, λ₀)

    @test pdae == pdae1
    @test pdae == pdae2

    @test hash(pdae1) == hash(pdae2)

    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)


    idae  = IDAE(eltype(q₀), 1, 1, 1, 1, pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, t₀, q₀, p₀, λ₀)
    idae1 = IDAE(pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, t₀, q₀, p₀, λ₀)
    idae2 = IDAE(pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, q₀, p₀, λ₀)

    @test idae == idae1
    @test idae == idae2

    @test hash(idae1) == hash(idae2)

    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)


    ################################################################################
    # Test HDAE: Hamiltonian Differential Algebraic Equation
    ################################################################################

    hdae_eqs = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_u, pdae_g, pdae_ϕ, pdae_ψ, pdae_h)

    hdae  = HDAE(eltype(q₀), 1, 1, 1, 1, hdae_eqs..., t₀, q₀, p₀, λ₀)
    hdae1 = HDAE(hdae_eqs..., t₀, q₀, p₀, λ₀)
    hdae2 = HDAE(hdae_eqs..., t₀, q₀, p₀)
    hdae3 = HDAE(hdae_eqs..., q₀, p₀)

    @test ndims(hdae) == 1
    @test periodicity(hdae) == zero(q₀)
    @test get_function_tuple(hdae) == NamedTuple{(:v, :f, :u, :g, :u̅, :g̅, :ϕ, :ψ, :h)}(hdae_eqs)

    @test hdae == hdae1
    @test hdae == hdae2
    @test hdae == hdae3

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)
    @test hash(hdae) == hash(hdae3)

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)


    ################################################################################
    # Test VDAE: Variational Differential Algebraic Equation
    ################################################################################

    function vdae_ψ(t, q, p, λ, μ, ψ)
        ψ[1] = μ[1] - λ[1]
    end

    vdae  = VDAE(eltype(q₀), 1, 1, 1, 1, iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, vdae_ψ, t₀, q₀, p₀, λ₀, λ₀; v=iode_v)
    vdae1 = VDAE(iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, vdae_ψ, t₀, q₀, p₀, λ₀, λ₀; v=iode_v)
    vdae2 = VDAE(iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, vdae_ψ, t₀, q₀, p₀, λ₀; v=iode_v)
    vdae3 = VDAE(iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, vdae_ψ, t₀, q₀, p₀; v=iode_v)
    vdae4 = VDAE(iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, vdae_ψ, q₀, p₀; v=iode_v)

    @test vdae == vdae1
    @test vdae == vdae2
    @test vdae == vdae3
    @test vdae == vdae4

    @test hash(vdae1) == hash(vdae2)
    @test hash(vdae3) == hash(vdae4)

    @test vdae == similar(vdae, t₀, q₀, p₀, λ₀, λ₀)
    @test vdae == similar(vdae, t₀, q₀, p₀, λ₀)
    @test vdae == similar(vdae, t₀, q₀, p₀)
    @test vdae == similar(vdae, q₀, p₀)


    ################################################################################
    # Test SODE: Split Ordinary Differential Equation
    ################################################################################

    function f_sode_1(t, x, f)
        f[1] = x[1]
    end

    function f_sode_2(t, x, f)
        f[1] = x[1]^2
    end

    f_sode = (f_sode_1, f_sode_2)

    sode  = SODE(eltype(q₀), 1, 1, 1, f_sode, t₀, q₀)
    sode1 = SODE(f_sode, t₀, q₀)
    sode2 = SODE(f_sode, q₀)

    @test sode == sode1
    @test sode == sode2

    @test hash(sode1) == hash(sode2)

    @test sode == similar(sode, t₀, q₀)
    @test sode == similar(sode, q₀)

end


using GeometricIntegrators.Common
using GeometricIntegrators.Equations
using Test

include("initial_conditions.jl")


function ode_v(t, x, ẋ)
    ẋ[1] = x[2]
    ẋ[2] = 2x[1]
end


function v_sode_1(t, x, v)
    v[1] = x[2]
end

function v_sode_2(t, x, v)
    v[2] = 2x[1]
end

function q_sode_1(t, x̄, x)
    x[1] = x̄[1]
    x[2] = x̄[2]
end

function q_sode_2(t, x̄, x)
    x[1] = x̄[1]
    x[2] = x̄[2]
end


function pode_v(t, q, p, v)
    v[1] = p[1]
end

function pode_f(t, q, p, f)
    f[1] = 2q[1]
end

function pode_h(t, q, p)
    p[1]^2/2 + cos(q[1])
end


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

function iode_h(t, q, v)
    v[1]^2/2 + cos(q[1])
end


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

function pdae_ψ(t, q, p, λ, μ, ψ)
    ψ[1] = μ[1] - λ[1]
end

function pdae_h(t, q, p)
    p[1]^2/2 + q[1]^2/2
end



@testset "$(rpad("Ordinary Differential Equations (ODE)",80))" begin

    ode  = ODE(ode_v, t₀, [x₀])
    ode1 = ODE(ode_v, [x₀])
    ode2 = ODE(ode_v, t₀, x₀)
    ode3 = ODE(ode_v, x₀)

    @test ndims(ode) == 2
    @test nsamples(ode) == 1
    @test periodicity(ode) == zero(x₀)
    @test get_function_tuple(ode) == NamedTuple{(:v,)}((ode_v,))

    @test ode == ode1
    @test ode == ode2
    @test ode == ode3

    @test hash(ode) == hash(ode1)
    @test hash(ode) == hash(ode2)
    @test hash(ode) == hash(ode3)

    @test ode == similar(ode, t₀, x₀)
    @test ode == similar(ode, x₀)

end


@testset "$(rpad("Split Ordinary Differential Equations (SODE)",80))" begin

    v_sode = (v_sode_1, v_sode_2)
    q_sode = (q_sode_1, q_sode_2)

    sode  = SODE(v_sode, t₀, [x₀])
    sode1 = SODE(v_sode, [x₀])
    sode2 = SODE(v_sode, t₀, x₀)
    sode3 = SODE(v_sode, x₀)

    @test ndims(sode) == 2
    @test nsamples(sode) == 1
    @test periodicity(sode) == zero(x₀)

    @test sode == sode1
    @test sode == sode2
    @test sode == sode3

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)

    @test sode == similar(sode, t₀, x₀)
    @test sode == similar(sode, x₀)


    sode  = SODE(v_sode, q_sode, t₀, [x₀])
    sode1 = SODE(v_sode, q_sode, [x₀])
    sode2 = SODE(v_sode, q_sode, t₀, x₀)
    sode3 = SODE(v_sode, q_sode, x₀)

    @test sode == sode1
    @test sode == sode2
    @test sode == sode3

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)

    @test sode == similar(sode, t₀, x₀)
    @test sode == similar(sode, x₀)

end


@testset "$(rpad("Partitioned Ordinary Differential Equations (PODE)",80))" begin

    pode_eqs = (pode_v, pode_f)

    pode  = PODE(pode_eqs..., t₀, [q₀], [p₀])
    pode1 = PODE(pode_eqs..., [q₀], [p₀])
    pode2 = PODE(pode_eqs..., t₀, q₀, p₀)
    pode3 = PODE(pode_eqs..., q₀, p₀)

    @test ndims(pode) == 1
    @test nsamples(pode) == 1
    @test periodicity(pode) == zero(q₀)
    @test get_function_tuple(pode) == NamedTuple{(:v, :f)}(pode_eqs)

    @test pode == pode1
    @test pode == pode2
    @test pode == pode3

    @test hash(pode) == hash(pode1)
    @test hash(pode) == hash(pode2)
    @test hash(pode) == hash(pode3)

    @test pode == similar(pode, t₀, q₀, p₀)
    @test pode == similar(pode, q₀, p₀)

    rode = ODE(ode_v, t₀, [x₀])
    code = convert(ODE, pode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v(rode.t₀, rode.q₀[begin], v₁)
    code.v(code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

    rode = SODE((v_sode_1, v_sode_2), t₀, [x₀])
    code = convert(SODE, pode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v[1](rode.t₀, rode.q₀[begin], v₁)
    code.v[1](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂
    rode.v[2](rode.t₀, rode.q₀[begin], v₁)
    code.v[2](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

end


@testset "$(rpad("Implicit Ordinary Differential Equations (IODE)",80))" begin

    iode_eqs = (iode_ϑ, iode_f, iode_g)

    # test without Hamiltonian

    iode  = IODE(iode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    iode1 = IODE(iode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    iode2 = IODE(iode_eqs..., t₀, q₀, p₀; v̄=iode_v)
    iode3 = IODE(iode_eqs..., q₀, p₀; v̄=iode_v)

    @test ndims(iode) == 1
    @test nsamples(iode) == 1
    @test periodicity(iode) == zero(q₀)
    @test get_function_tuple(iode) == NamedTuple{(:ϑ, :f, :g, :v̄, :f̄)}((iode_eqs..., iode_v, iode_f))

    @test iode == iode1
    @test iode == iode2
    @test iode == iode3

    @test hash(iode) == hash(iode1)
    @test hash(iode) == hash(iode2)
    @test hash(iode) == hash(iode3)

    @test iode == similar(iode, t₀, q₀, p₀, λ₀)
    @test iode == similar(iode, t₀, q₀, p₀)
    @test iode == similar(iode, q₀, p₀)
    

    # test with Hamiltonian

    iode  = IODE(iode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v, h=iode_h)
    iode1 = IODE(iode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v, h=iode_h)
    iode2 = IODE(iode_eqs..., t₀, q₀, p₀; v̄=iode_v, h=iode_h)
    iode3 = IODE(iode_eqs..., q₀, p₀; v̄=iode_v, h=iode_h)

    @test ndims(iode) == 1
    @test nsamples(iode) == 1
    @test periodicity(iode) == zero(q₀)
    @test get_function_tuple(iode) == NamedTuple{(:ϑ, :f, :g, :v̄, :f̄, :h)}((iode_eqs..., iode_v, iode_f, iode_h))

    @test iode == iode1
    @test iode == iode2
    @test iode == iode3

    @test hash(iode1) == hash(iode2)

    @test iode == similar(iode, t₀, q₀, p₀, λ₀)
    @test iode == similar(iode, t₀, q₀, p₀)
    @test iode == similar(iode, q₀, p₀)

end


@testset "$(rpad("Hamiltonian Ordinary Differential Equations (HODE)",80))" begin

    hode_eqs = (pode_v, pode_f, pode_h)

    hode  = HODE(hode_eqs..., t₀, [q₀], [p₀])
    hode1  = HODE(hode_eqs..., [q₀], [p₀])
    hode2 = HODE(hode_eqs..., t₀, q₀, p₀)
    hode3 = HODE(hode_eqs..., q₀, p₀)

    @test ndims(hode) == 1
    @test nsamples(hode) == 1
    @test periodicity(hode) == zero(q₀)
    @test get_function_tuple(hode) == NamedTuple{(:v, :f, :h)}(hode_eqs)

    @test hode == hode1
    @test hode == hode2
    @test hode == hode3

    @test hash(hode) == hash(hode1)
    @test hash(hode) == hash(hode2)
    @test hash(hode) == hash(hode3)

    @test hode == similar(hode, t₀, q₀, p₀)
    @test hode == similar(hode, q₀, p₀)

    rode = ODE(ode_v, t₀, [x₀])
    code = convert(ODE, hode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v(rode.t₀, rode.q₀[begin], v₁)
    code.v(code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

    pode = convert(PODE, hode)
    v₁ = zero(q₀)
    v₂ = zero(q₀)
    f₁ = zero(p₀)
    f₂ = zero(p₀)
    hode.v(hode.t₀, hode.q₀[begin], hode.p₀[begin], v₁)
    pode.v(pode.t₀, pode.q₀[begin], pode.p₀[begin], v₂)
    hode.f(hode.t₀, hode.q₀[begin], hode.p₀[begin], f₁)
    pode.f(pode.t₀, pode.q₀[begin], pode.p₀[begin], f₂)
    @test v₁ == v₂
    @test f₁ == f₂

    rode = SODE((v_sode_1, v_sode_2), t₀, [x₀])
    code = convert(SODE, hode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v[1](rode.t₀, rode.q₀[begin], v₁)
    code.v[1](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂
    rode.v[2](rode.t₀, rode.q₀[begin], v₁)
    code.v[2](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

end


@testset "$(rpad("Variational Ordinary Differential Equations (LODE)",80))" begin

    lode_eqs = (iode_ϑ, iode_f, iode_g)

    lode  = LODE(lode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    lode1 = LODE(lode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    lode2 = LODE(lode_eqs..., t₀, q₀, p₀; v̄=iode_v)
    lode3 = LODE(lode_eqs..., q₀, p₀; v̄=iode_v)

    @test ndims(lode) == 1
    @test nsamples(lode) == 1
    @test periodicity(lode) == zero(q₀)
    @test get_function_tuple(lode) == NamedTuple{(:ϑ, :f, :g, :v̄, :f̄)}((lode_eqs..., iode_v, iode_f))

    @test lode == lode1
    @test lode == lode2
    @test lode == lode3

    @test hash(lode) == hash(lode1)
    @test hash(lode) == hash(lode2)
    @test hash(lode) == hash(lode3)

    @test lode == similar(lode, t₀, q₀, p₀, λ₀)
    @test lode == similar(lode, t₀, q₀, p₀)
    @test lode == similar(lode, q₀, p₀)

end


@testset "$(rpad("Differential Algebraic Equations (DAE)",80))" begin

    dae_eqs = (dae_v, dae_u, dae_ϕ)

    dae  = DAE(dae_eqs..., t₀, [x₀], [λ₀])
    dae1 = DAE(dae_eqs..., [x₀], [λ₀])
    dae2 = DAE(dae_eqs..., t₀, x₀, λ₀)
    dae3 = DAE(dae_eqs..., x₀, λ₀)

    @test ndims(dae) == 2
    @test nsamples(dae) == 1
    @test nconstraints(dae) == 1
    @test periodicity(dae) == zero(x₀)
    @test get_function_tuple(dae) == NamedTuple{(:v, :u, :ϕ, :v̄)}((dae_eqs..., dae_v))

    @test dae == dae1
    @test dae == dae2
    @test dae == dae3

    @test hash(dae) == hash(dae1)
    @test hash(dae) == hash(dae2)
    @test hash(dae) == hash(dae3)

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)

end


@testset "$(rpad("Partitioned Differential Algebraic Equations (PDAE)",80))" begin

    pdae_eqs = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    pdae  = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    pdae1 = PDAE(pdae_eqs..., [q₀], [p₀], [λ₀])
    pdae2 = PDAE(pdae_eqs..., t₀, q₀, p₀, λ₀)
    pdae3 = PDAE(pdae_eqs..., q₀, p₀, λ₀)

    @test ndims(pdae) == 1
    @test nsamples(pdae) == 1
    @test nconstraints(pdae) == 1
    @test periodicity(pdae) == zero(q₀)
    @test get_function_tuple(pdae) == NamedTuple{(:v, :f, :u, :g, :ϕ, :v̄, :f̄)}((pdae_eqs..., pdae_v, pdae_f))

    @test pdae == pdae1
    @test pdae == pdae2
    @test pdae == pdae3

    @test hash(pdae) == hash(pdae1)
    @test hash(pdae) == hash(pdae2)
    @test hash(pdae) == hash(pdae3)

    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)

end


@testset "$(rpad("Implicit Differential Algebraic Equations (IDAE)",80))" begin

    idae_eqs = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    idae  = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae1 = IDAE(idae_eqs..., [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae2 = IDAE(idae_eqs..., t₀, q₀, p₀, λ₀; v̄=pdae_v)
    idae3 = IDAE(idae_eqs..., q₀, p₀, λ₀; v̄=pdae_v)

    @test ndims(idae) == 1
    @test nsamples(idae) == 1
    @test nconstraints(idae) == 1
    @test periodicity(idae) == zero(q₀)
    @test get_function_tuple(idae) == NamedTuple{(:ϑ, :f, :u, :g, :ϕ, :v̄, :f̄)}((idae_eqs..., pdae_v, pdae_f))

    @test idae == idae1
    @test idae == idae2
    @test idae == idae3

    @test hash(idae) == hash(idae1)
    @test hash(idae) == hash(idae2)
    @test hash(idae) == hash(idae3)

    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)

end


@testset "$(rpad("Hamiltonian Differential Algebraic Equations (HDAE)",80))" begin

    hdae_eqs = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_u, pdae_g, pdae_ϕ, pdae_ψ, pdae_h)

    hdae  = HDAE(hdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    hdae1 = HDAE(hdae_eqs..., [q₀], [p₀], [λ₀])
    hdae2 = HDAE(hdae_eqs..., t₀, q₀, p₀, λ₀)
    hdae3 = HDAE(hdae_eqs..., q₀, p₀, λ₀)

    @test ndims(hdae) == 1
    @test nsamples(hdae) == 1
    @test nconstraints(hdae) == 1
    @test periodicity(hdae) == zero(q₀)
    @test get_function_tuple(hdae) == NamedTuple{(:v, :f, :u, :g, :ū, :ḡ, :ϕ, :ψ, :h, :v̄, :f̄)}((hdae_eqs..., pdae_v, pdae_f))

    @test hdae == hdae1
    @test hdae == hdae2
    @test hdae == hdae3

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)
    @test hash(hdae) == hash(hdae3)

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)

end


@testset "$(rpad("Variational Differential Algebraic Equations (LDAE)",80))" begin

    ldae_eqs = (iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, pdae_ψ)

    ldae  = LDAE(ldae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
    ldae1 = LDAE(ldae_eqs..., [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
    ldae2 = LDAE(ldae_eqs..., t₀, q₀, p₀, λ₀, λ₀; v̄=iode_v)
    ldae3 = LDAE(ldae_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    ldae4 = LDAE(ldae_eqs..., t₀, q₀, p₀; v̄=iode_v)
    ldae5 = LDAE(ldae_eqs..., q₀, p₀; v̄=iode_v)

    @test ndims(ldae) == 1
    @test nsamples(ldae) == 1
    @test nconstraints(ldae) == 1
    @test periodicity(ldae) == zero(q₀)
    @test get_function_tuple(ldae) == NamedTuple{(:ϑ, :f, :g, :ḡ, :ϕ, :ψ, :v̄, :f̄)}((ldae_eqs..., iode_v, iode_f))

    @test ldae == ldae1
    @test ldae == ldae2
    @test ldae == ldae3
    @test ldae == ldae4
    @test ldae == ldae5

    @test hash(ldae) == hash(ldae1)
    @test hash(ldae) == hash(ldae2)
    @test hash(ldae) == hash(ldae3)
    @test hash(ldae) == hash(ldae4)
    @test hash(ldae) == hash(ldae5)

    @test ldae == similar(ldae, t₀, q₀, p₀, λ₀, λ₀)
    @test ldae == similar(ldae, t₀, q₀, p₀, λ₀)
    @test ldae == similar(ldae, t₀, q₀, p₀)
    @test ldae == similar(ldae, q₀, p₀)

end


@testset "$(rpad("Split Partitioned Differential Algebraic Equations (SPDAE)",80))" begin

    spdae_eqs = ((pdae_v, pdae_u, pdae_u), (pdae_f, pdae_g, pdae_g), pdae_ϕ, pdae_ψ)

    spdae  = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    spdae1 = SPDAE(spdae_eqs..., [q₀], [p₀], [λ₀])
    spdae2 = SPDAE(spdae_eqs..., t₀, q₀, p₀, λ₀)
    spdae3 = SPDAE(spdae_eqs..., t₀, q₀, p₀)
    spdae4 = SPDAE(spdae_eqs..., q₀, p₀)

    @test ndims(spdae) == 1
    @test nsamples(spdae) == 1
    @test nconstraints(spdae) == 1
    @test periodicity(spdae) == zero(q₀)
    @test get_function_tuple(spdae) == NamedTuple{(:v, :f, :ϕ, :ψ)}(spdae_eqs)

    @test spdae == spdae1
    @test spdae == spdae2
    @test spdae == spdae3
    @test spdae == spdae4

    @test hash(spdae) == hash(spdae1)
    @test hash(spdae) == hash(spdae2)
    @test hash(spdae) == hash(spdae3)
    @test hash(spdae) == hash(spdae4)

    @test spdae == similar(spdae, t₀, q₀, p₀, λ₀)
    @test spdae == similar(spdae, t₀, q₀, p₀)
    @test spdae == similar(spdae, q₀, p₀)

end

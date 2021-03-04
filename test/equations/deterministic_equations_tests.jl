
using GeometricIntegrators.Common
using GeometricIntegrators.Equations
using GeometricIntegrators.Utils
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

function iode_u(t, q, λ, u)
    u[1] = λ[1]
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

function lode_l(t, q, v)
    v[1]^2/2 - cos(q[1])
end

function lode_ω(t, q, v, ω)
    ω[1,1] = sin(q[1])
    ω[1,2] = 0
    ω[2,1] = 0
    ω[2,2] = v[1]
end


function dae_v(t, x, v)
    v[1] = x[1]
    v[2] = x[2]
end

function dae_u(t, x, λ, u)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ū(t, x, λ, u)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ϕ(t, x, λ, ϕ)
    ϕ[1] = x[2] - x[1]
end

function dae_ψ(t, x, v, λ, ψ)
    ψ[1] = v[2] - v[1]
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

    ode  = ODE(ode_v, t₀, [x₀], nothing, nothing, nothing)

    ode1 = ODE(ode_v, t₀, [x₀])
    ode2 = ODE(ode_v, [x₀])
    ode3 = ODE(ode_v, t₀, x₀)
    ode4 = ODE(ode_v, x₀)

    @test axes(ode) == axes(x₀)
    @test ndims(ode) == 2
    @test nsamples(ode) == 1

    @test periodicity(ode) == zero(x₀)
    @test initial_conditions(ode) == (t₀, [x₀])

    @test hasinvariants(ode) == false
    @test hasparameters(ode) == false
    @test hasperiodicity(ode) == false

    @test get_function_tuple(ode) == NamedTuple{(:v,)}((ode_v,))

    @test ode == ode1
    @test ode == ode2
    @test ode == ode3
    @test ode == ode4

    @test hash(ode) == hash(ode1)
    @test hash(ode) == hash(ode2)
    @test hash(ode) == hash(ode3)
    @test hash(ode) == hash(ode4)

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

    pode  = PODE(pode_eqs..., t₀, [q₀], [p₀], nothing, nothing, nothing)

    pode1 = PODE(pode_eqs..., t₀, [q₀], [p₀])
    pode2 = PODE(pode_eqs..., [q₀], [p₀])
    pode3 = PODE(pode_eqs..., t₀, q₀, p₀)
    pode4 = PODE(pode_eqs..., q₀, p₀)

    @test axes(pode) == axes(q₀)
    @test ndims(pode) == 1
    @test nsamples(pode) == 1

    @test periodicity(pode) == zero(q₀)
    @test initial_conditions(pode) == (t₀, [q₀], [p₀])

    @test hasinvariants(pode) == false
    @test hasparameters(pode) == false
    @test hasperiodicity(pode) == false

    functions = get_function_tuple(pode)
    @test functions.v == pode_v == pode.v
    @test functions.f == pode_f == pode.f

    @test pode == pode1
    @test pode == pode2
    @test pode == pode3
    @test pode == pode4

    @test hash(pode) == hash(pode1)
    @test hash(pode) == hash(pode2)
    @test hash(pode) == hash(pode3)
    @test hash(pode) == hash(pode4)

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

    @test axes(iode) == axes(q₀)
    @test ndims(iode) == 1
    @test nsamples(iode) == 1

    @test periodicity(iode) == zero(q₀)
    @test initial_conditions(iode) == (t₀, [q₀], [p₀], [λ₀])

    @test hasinvariants(iode) == false
    @test hasparameters(iode) == false
    @test hasperiodicity(iode) == false

    functions = get_function_tuple(iode)
    @test functions.ϑ == iode_ϑ == iode.ϑ
    @test functions.f == iode_f == iode.f
    @test functions.g == iode_g == iode.g
    @test functions.v̄ == iode_v == iode.v̄
    @test functions.f̄ == iode_f == iode.f̄

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

    # iode  = IODE(iode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v, h=iode_h)
    # iode1 = IODE(iode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v, h=iode_h)
    # iode2 = IODE(iode_eqs..., t₀, q₀, p₀; v̄=iode_v, h=iode_h)
    # iode3 = IODE(iode_eqs..., q₀, p₀; v̄=iode_v, h=iode_h)

    # @test ndims(iode) == 1
    # @test nsamples(iode) == 1
    # @test periodicity(iode) == zero(q₀)
    # @test get_function_tuple(iode) == NamedTuple{(:ϑ, :f, :g, :v̄, :f̄, :h)}((iode_eqs..., iode_v, iode_f, iode_h))

    # @test iode == iode1
    # @test iode == iode2
    # @test iode == iode3

    # @test hash(iode1) == hash(iode2)

    # @test iode == similar(iode, t₀, q₀, p₀, λ₀)
    # @test iode == similar(iode, t₀, q₀, p₀)
    # @test iode == similar(iode, q₀, p₀)

end


@testset "$(rpad("Hamiltonian Ordinary Differential Equations (HODE)",80))" begin

    hode_eqs = (pode_v, pode_f, pode_h)

    hode  = HODE(hode_eqs..., t₀, [q₀], [p₀])
    hode1  = HODE(hode_eqs..., [q₀], [p₀])
    hode2 = HODE(hode_eqs..., t₀, q₀, p₀)
    hode3 = HODE(hode_eqs..., q₀, p₀)

    @test axes(hode) == axes(q₀)
    @test ndims(hode) == 1
    @test nsamples(hode) == 1

    @test periodicity(hode) == zero(q₀)
    @test initial_conditions(hode) == (t₀, [q₀], [p₀])

    @test hasinvariants(hode) == false
    @test hasparameters(hode) == false
    @test hasperiodicity(hode) == false

    functions = get_function_tuple(hode)
    @test functions.v == pode_v == hode.v
    @test functions.f == pode_f == hode.f
    @test functions.h == pode_h == hode.hamiltonian

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

    lode_eqs = (iode_ϑ, iode_f, iode_g, lode_l, lode_ω)

    lode  = LODE(iode_ϑ, iode_f, iode_g, lode_ω, iode_v, iode_f, t₀, [q₀], [p₀], [λ₀], lode_l, nothing, nothing, nothing)

    lode1 = LODE(lode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    lode2 = LODE(lode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    lode3 = LODE(lode_eqs..., t₀, q₀, p₀; v̄=iode_v)
    lode4 = LODE(lode_eqs..., q₀, p₀; v̄=iode_v)

    @test axes(lode) == axes(q₀)
    @test ndims(lode) == 1
    @test nsamples(lode) == 1

    @test periodicity(lode) == zero(q₀)
    @test initial_conditions(lode) == (t₀, [q₀], [p₀], [λ₀])

    @test hasinvariants(lode) == false
    @test hasparameters(lode) == false
    @test hasperiodicity(lode) == false

    functions = get_function_tuple(lode)
    @test functions.ϑ == iode_ϑ == lode.ϑ
    @test functions.f == iode_f == lode.f
    @test functions.g == iode_g == lode.g
    @test functions.ω == lode_ω == lode.ω
    @test functions.v̄ == iode_v == lode.v̄
    @test functions.f̄ == iode_f == lode.f̄
    @test functions.l == lode_l == lode.lagrangian

    @test lode == lode1
    @test lode == lode2
    @test lode == lode3
    @test lode == lode4

    @test hash(lode) == hash(lode1)
    @test hash(lode) == hash(lode2)
    @test hash(lode) == hash(lode3)
    @test hash(lode) == hash(lode4)

    @test lode == similar(lode, t₀, q₀, p₀, λ₀)
    @test lode == similar(lode, t₀, q₀, p₀)
    @test lode == similar(lode, q₀, p₀)

    iode = convert(IODE, lode)
    v₁ = zero(q₀)
    v₂ = zero(q₀)
    ϑ₁ = zero(q₀)
    ϑ₂ = zero(q₀)
    f₁ = zero(p₀)
    f₂ = zero(p₀)
    g₁ = zero(p₀)
    g₂ = zero(p₀)
    lode.v̄(lode.t₀, lode.q₀[begin], lode.p₀[begin], v₁)
    iode.v̄(iode.t₀, iode.q₀[begin], iode.p₀[begin], v₂)
    lode.ϑ(lode.t₀, lode.q₀[begin], v₁, ϑ₁)
    iode.ϑ(iode.t₀, iode.q₀[begin], v₂, ϑ₂)
    lode.f(lode.t₀, lode.q₀[begin], v₁, f₁)
    iode.f(iode.t₀, iode.q₀[begin], v₂, f₂)
    lode.g(lode.t₀, lode.q₀[begin], v₁, g₁)
    iode.g(iode.t₀, iode.q₀[begin], v₂, g₂)
    @test v₁ == v₂
    @test ϑ₁ == ϑ₂
    @test f₁ == f₂
    @test g₁ == g₂

end


@testset "$(rpad("Differential Algebraic Equations (DAE)",80))" begin

    dae_eqs  = (dae_v, dae_u, nothing, dae_ϕ, nothing, dae_v)
    dae_eqs1 = (dae_v, dae_u, dae_ϕ)
    dae_eqs2 = (dae_v, dae_u, nothing, dae_ϕ, nothing)

    dae  = DAE(dae_eqs..., t₀, [x₀], [λ₀], nothing, nothing, nothing)

    dae1 = DAE(dae_eqs1..., [x₀], [λ₀])
    dae2 = DAE(dae_eqs1..., [x₀], [λ₀])
    dae3 = DAE(dae_eqs1..., t₀, x₀, λ₀)
    dae4 = DAE(dae_eqs1..., x₀, λ₀)

    dae5 = DAE(dae_eqs2..., [x₀], [λ₀])
    dae6 = DAE(dae_eqs2..., [x₀], [λ₀])
    dae7 = DAE(dae_eqs2..., t₀, x₀, λ₀)
    dae8 = DAE(dae_eqs2..., x₀, λ₀)
    
    @test axes(dae) == axes(x₀)
    @test ndims(dae) == 2
    @test nsamples(dae) == 1
    @test nconstraints(dae) == 1

    @test periodicity(dae) == zero(x₀)
    @test initial_conditions(dae) == (t₀, [x₀], [λ₀])

    @test hassecondary(dae) == false
    @test hasinvariants(dae) == false
    @test hasparameters(dae) == false
    @test hasperiodicity(dae) == false

    functions = get_function_tuple(dae)
    @test functions.v == dae_v == dae.v
    @test functions.u == dae_u == dae.u
    @test functions.ϕ == dae_ϕ == dae.ϕ
    @test functions.v̄ == dae_v == dae.v̄

    @test dae == dae1
    @test dae == dae2
    @test dae == dae3
    @test dae == dae4
    @test dae == dae5
    @test dae == dae6
    @test dae == dae7
    @test dae == dae8

    @test hash(dae) == hash(dae1)
    @test hash(dae) == hash(dae2)
    @test hash(dae) == hash(dae3)
    @test hash(dae) == hash(dae4)
    @test hash(dae) == hash(dae5)
    @test hash(dae) == hash(dae6)
    @test hash(dae) == hash(dae7)
    @test hash(dae) == hash(dae8)

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)


    dae_eqs  = (dae_v, dae_u, dae_ū, dae_ϕ, dae_ψ)
    dae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(2))

    dae  = DAE(dae_eqs..., dae_v, t₀, [x₀], [λ₀], dae_args.invariants, dae_args.parameters, dae_args.periodicity)

    dae1 = DAE(dae_eqs..., [x₀], [λ₀]; dae_args...)
    dae2 = DAE(dae_eqs..., [x₀], [λ₀]; dae_args...)
    dae3 = DAE(dae_eqs..., t₀, x₀, λ₀; dae_args...)
    dae4 = DAE(dae_eqs..., x₀, λ₀; dae_args...)

    @test axes(dae) == axes(x₀)
    @test ndims(dae) == 2
    @test nsamples(dae) == 1
    @test nconstraints(dae) == 1

    @test periodicity(dae) == dae_args.periodicity
    @test initial_conditions(dae) == (t₀, [x₀], [λ₀])

    @test hassecondary(dae) == true
    @test hasinvariants(dae) == true
    @test hasparameters(dae) == true
    @test hasperiodicity(dae) == true

    functions = get_function_tuple(dae)
    @test functions.v != dae_v == dae.v
    @test functions.u != dae_u == dae.u
    @test functions.ū != dae_ū == dae.ū
    @test functions.ϕ != dae_ϕ == dae.ϕ
    @test functions.ψ != dae_ψ == dae.ψ
    @test functions.v̄ != dae_v == dae.v̄

    @test dae == dae1
    @test dae == dae2
    @test dae == dae3
    @test dae == dae4

    @test hash(dae) == hash(dae1)
    @test hash(dae) == hash(dae2)
    @test hash(dae) == hash(dae3)
    @test hash(dae) == hash(dae4)

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)

end


@testset "$(rpad("Partitioned Differential Algebraic Equations (PDAE)",80))" begin

    pdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f)
    pdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)
    pdae_eqs2 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing)

    pdae  = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀], nothing, nothing, nothing)

    pdae1 = PDAE(pdae_eqs1..., t₀, [q₀], [p₀], [λ₀])
    pdae2 = PDAE(pdae_eqs1..., [q₀], [p₀], [λ₀])
    pdae3 = PDAE(pdae_eqs1..., t₀, q₀, p₀, λ₀)
    pdae4 = PDAE(pdae_eqs1..., q₀, p₀, λ₀)

    pdae5 = PDAE(pdae_eqs2..., t₀, [q₀], [p₀], [λ₀])
    pdae6 = PDAE(pdae_eqs2..., [q₀], [p₀], [λ₀])
    pdae7 = PDAE(pdae_eqs2..., t₀, q₀, p₀, λ₀)
    pdae8 = PDAE(pdae_eqs2..., q₀, p₀, λ₀)

    @test axes(pdae) == axes(q₀)
    @test ndims(pdae) == 1
    @test nsamples(pdae) == 1
    @test nconstraints(pdae) == 1

    @test periodicity(pdae) == zero(q₀)
    @test initial_conditions(pdae) == (t₀, [q₀], [p₀], [λ₀])

    @test hassecondary(pdae) == false
    @test hasinvariants(pdae) == false
    @test hasparameters(pdae) == false
    @test hasperiodicity(pdae) == false

    functions = get_function_tuple(pdae)
    @test functions.v == pdae_v == pdae.v
    @test functions.f == pdae_f == pdae.f
    @test functions.u == pdae_u == pdae.u
    @test functions.g == pdae_g == pdae.g
    @test functions.ϕ == pdae_ϕ == pdae.ϕ
    @test functions.v̄ == pdae_v == pdae.v̄
    @test functions.f̄ == pdae_f == pdae.f̄
    
    @test pdae == pdae1
    @test pdae == pdae2
    @test pdae == pdae3
    @test pdae == pdae4
    @test pdae == pdae5
    @test pdae == pdae6
    @test pdae == pdae7
    @test pdae == pdae8

    @test hash(pdae) == hash(pdae1)
    @test hash(pdae) == hash(pdae2)
    @test hash(pdae) == hash(pdae3)
    @test hash(pdae) == hash(pdae4)
    @test hash(pdae) == hash(pdae5)
    @test hash(pdae) == hash(pdae6)
    @test hash(pdae) == hash(pdae7)
    @test hash(pdae) == hash(pdae8)

    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)


    pdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    pdae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    pdae  = PDAE(pdae_eqs..., pdae_v, pdae_f, t₀, [q₀], [p₀], [λ₀], pdae_args.invariants, pdae_args.parameters, pdae_args.periodicity)

    pdae1 = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀]; pdae_args...)
    pdae2 = PDAE(pdae_eqs..., [q₀], [p₀], [λ₀]; pdae_args...)
    pdae3 = PDAE(pdae_eqs..., t₀, q₀, p₀, λ₀; pdae_args...)
    pdae4 = PDAE(pdae_eqs..., q₀, p₀, λ₀; pdae_args...)

    @test axes(pdae) == axes(q₀)
    @test ndims(pdae) == 1
    @test nsamples(pdae) == 1
    @test nconstraints(pdae) == 1

    @test periodicity(pdae) == pdae_args.periodicity
    @test initial_conditions(pdae) == (t₀, [q₀], [p₀], [λ₀])

    @test hassecondary(pdae) == true
    @test hasinvariants(pdae) == true
    @test hasparameters(pdae) == true
    @test hasperiodicity(pdae) == true

    functions = get_function_tuple(pdae)
    @test functions.v != pdae_v == pdae.v
    @test functions.f != pdae_f == pdae.f
    @test functions.u != pdae_u == pdae.u
    @test functions.g != pdae_g == pdae.g
    @test functions.ϕ != pdae_ϕ == pdae.ϕ
    @test functions.ū != pdae_u == pdae.ū
    @test functions.ḡ != pdae_g == pdae.ḡ
    @test functions.ψ != pdae_ψ == pdae.ψ
    @test functions.v̄ != pdae_v == pdae.v̄
    @test functions.f̄ != pdae_f == pdae.f̄
    
    @test pdae == pdae1
    @test pdae == pdae2
    @test pdae == pdae3
    @test pdae == pdae4

    @test hash(pdae) == hash(pdae1)
    @test hash(pdae) == hash(pdae2)
    @test hash(pdae) == hash(pdae3)
    @test hash(pdae) == hash(pdae4)

    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)

end


@testset "$(rpad("Implicit Differential Algebraic Equations (IDAE)",80))" begin

    idae_eqs  = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f)
    idae_eqs1 = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ)
    idae_eqs2 = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing)

    idae  = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀], nothing, nothing, nothing)

    idae1 = IDAE(idae_eqs1..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae2 = IDAE(idae_eqs1..., [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae3 = IDAE(idae_eqs1..., t₀, q₀, p₀, λ₀; v̄=pdae_v)
    idae4 = IDAE(idae_eqs1..., q₀, p₀, λ₀; v̄=pdae_v)

    idae5 = IDAE(idae_eqs2..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae6 = IDAE(idae_eqs2..., [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae7 = IDAE(idae_eqs2..., t₀, q₀, p₀, λ₀; v̄=pdae_v)
    idae8 = IDAE(idae_eqs2..., q₀, p₀, λ₀; v̄=pdae_v)

    @test axes(idae) == axes(q₀)
    @test ndims(idae) == 1
    @test nsamples(idae) == 1
    @test nconstraints(idae) == 1

    @test periodicity(idae) == zero(q₀)
    @test initial_conditions(idae) == (t₀, [q₀], [p₀], [λ₀])

    @test hassecondary(idae) == false
    @test hasinvariants(idae) == false
    @test hasparameters(idae) == false
    @test hasperiodicity(idae) == false

    functions = get_function_tuple(idae)
    @test functions.ϑ == pdae_p == idae.ϑ
    @test functions.f == pdae_f == idae.f
    @test functions.u == pdae_u == idae.u
    @test functions.g == pdae_g == idae.g
    @test functions.ϕ == pdae_ϕ == idae.ϕ
    @test functions.v̄ == pdae_v == idae.v̄
    @test functions.f̄ == pdae_f == idae.f̄
    
    @test idae == idae1
    @test idae == idae2
    @test idae == idae3
    @test idae == idae4
    @test idae == idae5
    @test idae == idae6
    @test idae == idae7
    @test idae == idae8

    @test hash(idae) == hash(idae1)
    @test hash(idae) == hash(idae2)
    @test hash(idae) == hash(idae3)
    @test hash(idae) == hash(idae4)
    @test hash(idae) == hash(idae5)
    @test hash(idae) == hash(idae6)
    @test hash(idae) == hash(idae7)
    @test hash(idae) == hash(idae8)

    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)


    idae_eqs  = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    idae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    idae  = IDAE(idae_eqs..., pdae_v, pdae_f, t₀, [q₀], [p₀], [λ₀], [λ₀], idae_args.invariants, idae_args.parameters, idae_args.periodicity)

    idae1 = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v, idae_args...)
    idae2 = IDAE(idae_eqs..., [q₀], [p₀], [λ₀]; v̄=pdae_v, idae_args...)
    idae3 = IDAE(idae_eqs..., t₀, q₀, p₀, λ₀; v̄=pdae_v, idae_args...)
    idae4 = IDAE(idae_eqs..., q₀, p₀, λ₀; v̄=pdae_v, idae_args...)

    @test axes(idae) == axes(q₀)
    @test ndims(idae) == 1
    @test nsamples(idae) == 1
    @test nconstraints(idae) == 1

    @test periodicity(idae) == idae_args.periodicity
    @test initial_conditions(idae) == (t₀, [q₀], [p₀], [λ₀])

    @test hassecondary(idae) == true
    @test hasinvariants(idae) == true
    @test hasparameters(idae) == true
    @test hasperiodicity(idae) == true

    functions = get_function_tuple(idae)
    @test functions.ϑ != pdae_p == idae.ϑ
    @test functions.f != pdae_f == idae.f
    @test functions.u != pdae_u == idae.u
    @test functions.g != pdae_g == idae.g
    @test functions.ϕ != pdae_ϕ == idae.ϕ
    @test functions.ū != pdae_u == idae.ū
    @test functions.ḡ != pdae_g == idae.ḡ
    @test functions.ψ != pdae_ψ == idae.ψ
    @test functions.v̄ != pdae_v == idae.v̄
    @test functions.f̄ != pdae_f == idae.f̄
    
    @test idae == idae1
    @test idae == idae2
    @test idae == idae3
    @test idae == idae4

    @test hash(idae) == hash(idae1)
    @test hash(idae) == hash(idae2)
    @test hash(idae) == hash(idae3)
    @test hash(idae) == hash(idae4)

    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)

end


@testset "$(rpad("Hamiltonian Differential Algebraic Equations (HDAE)",80))" begin

    hdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f, symplectic_matrix)
    hdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_h)
    hdae_eqs2 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_h)

    hdae  = HDAE(hdae_eqs..., t₀, [q₀], [p₀], [λ₀], pdae_h, nothing, nothing, nothing)

    hdae1 = HDAE(hdae_eqs1..., t₀, [q₀], [p₀], [λ₀])
    hdae2 = HDAE(hdae_eqs1..., [q₀], [p₀], [λ₀])
    hdae3 = HDAE(hdae_eqs1..., t₀, q₀, p₀, λ₀)
    hdae4 = HDAE(hdae_eqs1..., q₀, p₀, λ₀)

    hdae5 = HDAE(hdae_eqs2..., t₀, [q₀], [p₀], [λ₀])
    hdae6 = HDAE(hdae_eqs2..., [q₀], [p₀], [λ₀])
    hdae7 = HDAE(hdae_eqs2..., t₀, q₀, p₀, λ₀)
    hdae8 = HDAE(hdae_eqs2..., q₀, p₀, λ₀)

    @test axes(hdae) == axes(q₀)
    @test ndims(hdae) == 1
    @test nsamples(hdae) == 1
    @test nconstraints(hdae) == 1

    @test periodicity(hdae) == zero(q₀)
    @test initial_conditions(hdae) == (t₀, [q₀], [p₀], [λ₀])

    @test hassecondary(hdae) == false
    @test hasinvariants(hdae) == false
    @test hasparameters(hdae) == false
    @test hasperiodicity(hdae) == false

    functions = get_function_tuple(hdae)
    @test functions.v == pdae_v == hdae.v
    @test functions.f == pdae_f == hdae.f
    @test functions.u == pdae_u == hdae.u
    @test functions.g == pdae_g == hdae.g
    @test functions.ϕ == pdae_ϕ == hdae.ϕ
    @test functions.v̄ == pdae_v == hdae.v̄
    @test functions.f̄ == pdae_f == hdae.f̄
    @test functions.h == pdae_h == hdae.hamiltonian

    @test hdae == hdae1
    @test hdae == hdae2
    @test hdae == hdae3
    @test hdae == hdae4
    @test hdae == hdae5
    @test hdae == hdae6
    @test hdae == hdae7
    @test hdae == hdae8

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)
    @test hash(hdae) == hash(hdae3)
    @test hash(hdae) == hash(hdae4)
    @test hash(hdae) == hash(hdae5)
    @test hash(hdae) == hash(hdae6)
    @test hash(hdae) == hash(hdae7)
    @test hash(hdae) == hash(hdae8)

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)


    hdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    hdae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    hdae  = HDAE(hdae_eqs..., pdae_v, pdae_f, symplectic_matrix, t₀, [q₀], [p₀], [λ₀], pdae_h, hdae_args.invariants, hdae_args.parameters, hdae_args.periodicity)

    hdae1 = HDAE(hdae_eqs..., pdae_h, t₀, [q₀], [p₀], [λ₀]; hdae_args...)
    hdae2 = HDAE(hdae_eqs..., pdae_h, [q₀], [p₀], [λ₀]; hdae_args...)
    hdae3 = HDAE(hdae_eqs..., pdae_h, t₀, q₀, p₀, λ₀; hdae_args...)
    hdae4 = HDAE(hdae_eqs..., pdae_h, q₀, p₀, λ₀; hdae_args...)

    @test axes(hdae) == axes(q₀)
    @test ndims(hdae) == 1
    @test nsamples(hdae) == 1
    @test nconstraints(hdae) == 1

    @test periodicity(hdae) == hdae_args.periodicity
    @test initial_conditions(hdae) == (t₀, [q₀], [p₀], [λ₀])

    @test hassecondary(hdae) == true
    @test hasinvariants(hdae) == true
    @test hasparameters(hdae) == true
    @test hasperiodicity(hdae) == true

    functions = get_function_tuple(hdae)
    @test functions.v != pdae_v == hdae.v
    @test functions.f != pdae_f == hdae.f
    @test functions.u != pdae_u == hdae.u
    @test functions.g != pdae_g == hdae.g
    @test functions.ϕ != pdae_ϕ == hdae.ϕ
    @test functions.ū != pdae_u == hdae.ū
    @test functions.ḡ != pdae_g == hdae.ḡ
    @test functions.ψ != pdae_ψ == hdae.ψ
    @test functions.v̄ != pdae_v == hdae.v̄
    @test functions.f̄ != pdae_f == hdae.f̄
    @test functions.h != pdae_h == hdae.hamiltonian

    @test hdae == hdae1
    @test hdae == hdae2
    @test hdae == hdae3
    @test hdae == hdae4

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)
    @test hash(hdae) == hash(hdae3)
    @test hash(hdae) == hash(hdae4)

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)

end


@testset "$(rpad("Variational Differential Algebraic Equations (LDAE)",80))" begin

    ldae_eqs = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, nothing, nothing, nothing, lode_l, lode_ω)

    ldae  = LDAE(ldae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)

    ldae01 = LDAE(ldae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
    ldae02 = LDAE(ldae_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    ldae03 = LDAE(ldae_eqs..., [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
    ldae04 = LDAE(ldae_eqs..., [q₀], [p₀], [λ₀]; v̄=iode_v)
    ldae05 = LDAE(ldae_eqs..., t₀, q₀, p₀, λ₀, λ₀; v̄=iode_v)
    ldae06 = LDAE(ldae_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    ldae07 = LDAE(ldae_eqs..., q₀, p₀, λ₀, λ₀; v̄=iode_v)
    ldae08 = LDAE(ldae_eqs..., q₀, p₀, λ₀; v̄=iode_v)

    @test ndims(ldae) == 1
    @test nsamples(ldae) == 1
    @test nconstraints(ldae) == 1
    @test periodicity(ldae) == zero(q₀)

    @test periodicity(ldae) == zero(q₀)
    @test initial_conditions(ldae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(ldae) == false
    @test hasinvariants(ldae) == false
    @test hasparameters(ldae) == false
    @test hasperiodicity(ldae) == false

    functions = get_function_tuple(ldae)
    @test functions.ϑ == iode_ϑ == ldae.ϑ
    @test functions.f == iode_f == ldae.f
    @test functions.u == iode_u == ldae.u
    @test functions.g == iode_g == ldae.g
    @test functions.ϕ == pdae_ϕ == ldae.ϕ
    @test functions.l == lode_l == ldae.lagrangian
    @test functions.ω == lode_ω == ldae.ω
    @test functions.v̄ == iode_v == ldae.v̄
    @test functions.f̄ == iode_f == ldae.f̄
    
    @test ldae == ldae01
    @test ldae == ldae02
    @test ldae == ldae03
    @test ldae == ldae04
    @test ldae == ldae05
    @test ldae == ldae06
    @test ldae == ldae07
    @test ldae == ldae08
    # @test ldae == ldae09
    # @test ldae == ldae10
    # @test ldae == ldae11
    # @test ldae == ldae12

    @test hash(ldae) == hash(ldae01)
    @test hash(ldae) == hash(ldae02)
    @test hash(ldae) == hash(ldae03)
    @test hash(ldae) == hash(ldae04)
    @test hash(ldae) == hash(ldae05)
    @test hash(ldae) == hash(ldae06)
    @test hash(ldae) == hash(ldae07)
    @test hash(ldae) == hash(ldae08)
    # @test hash(ldae) == hash(ldae09)
    # @test hash(ldae) == hash(ldae10)
    # @test hash(ldae) == hash(ldae11)
    # @test hash(ldae) == hash(ldae12)

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

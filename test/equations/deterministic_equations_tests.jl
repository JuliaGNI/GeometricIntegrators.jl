
using GeometricIntegrators.Common
using GeometricIntegrators.Equations
using Test

include("initial_conditions.jl")


function ode_v(t, x, f)
    f[1] = x[1]
end


function pode_v(t, q, p, v)
    v[1] = q[1]
end

function pode_f(t, q, p, f)
    f[1] = 2p[1]
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


function f_sode_1(t, x, f)
    f[1] = x[1]
end

function f_sode_2(t, x, f)
    f[1] = x[1]^2
end

function q_sode_1(t, x̄, x)
    x[1] = x̄[1]
end

function q_sode_2(t, x̄, x)
    x[1] = x̄[1]
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

    ode  = ODE(ode_v, t₀, [q₀])
    ode1 = ODE(ode_v, t₀, q₀)
    ode2 = ODE(ode_v, q₀)

    @test ndims(ode) == 1
    @test nsamples(ode) == 1
    @test periodicity(ode) == zero(q₀)
    @test get_function_tuple(ode) == NamedTuple{(:v,)}((ode_v,))

    @test ode == ode1
    @test ode == ode2

    @test hash(ode1) == hash(ode2)

    @test ode == similar(ode, t₀, q₀)
    @test ode == similar(ode, q₀)

end


@testset "$(rpad("Partitioned Ordinary Differential Equations (PODE)",80))" begin

    pode_eqs = (pode_v, pode_f)

    pode  = PODE(pode_eqs..., t₀, [q₀], [p₀])
    pode1 = PODE(pode_eqs..., t₀, q₀, p₀)
    pode2 = PODE(pode_eqs..., q₀, p₀)

    @test ndims(pode) == 1
    @test nsamples(pode) == 1
    @test periodicity(pode) == zero(q₀)
    @test get_function_tuple(pode) == NamedTuple{(:v, :f)}(pode_eqs)

    @test pode == pode1
    @test pode == pode2

    @test hash(pode1) == hash(pode2)

    @test pode == similar(pode, t₀, q₀, p₀)
    @test pode == similar(pode, q₀, p₀)

end


@testset "$(rpad("Implicit Ordinary Differential Equations (IODE)",80))" begin

    iode_eqs = (iode_ϑ, iode_f, iode_g)

    # test without Hamiltonian

    iode  = IODE(eltype(q₀), 1, 1, 1, iode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
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

    @test hash(iode1) == hash(iode2)

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

    hode_eqs = (iode_v, iode_f, iode_h)

    hode  = HODE(hode_eqs..., t₀, [q₀], [p₀])
    hode1 = HODE(hode_eqs..., t₀, q₀, p₀)
    hode2 = HODE(hode_eqs..., q₀, p₀)

    @test ndims(hode) == 1
    @test nsamples(hode) == 1
    @test periodicity(hode) == zero(q₀)
    @test get_function_tuple(hode) == NamedTuple{(:v, :f, :h)}(hode_eqs)

    @test hode == hode1
    @test hode == hode2

    @test hash(hode1) == hash(hode2)

    @test hode == similar(hode, t₀, q₀, p₀)
    @test hode == similar(hode, q₀, p₀)

end


@testset "$(rpad("Variational Ordinary Differential Equations (VODE)",80))" begin

    vode_eqs = (iode_ϑ, iode_f, iode_g)

    vode  = VODE(vode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    vode1 = VODE(vode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    vode2 = VODE(vode_eqs..., t₀, q₀, p₀; v̄=iode_v)
    vode3 = VODE(vode_eqs..., q₀, p₀; v̄=iode_v)

    @test ndims(vode) == 1
    @test nsamples(vode) == 1
    @test periodicity(vode) == zero(q₀)
    @test get_function_tuple(vode) == NamedTuple{(:ϑ, :f, :g, :v̄, :f̄)}((vode_eqs..., iode_v, iode_f))

    @test vode == vode1
    @test vode == vode2
    @test vode == vode3

    @test hash(vode1) == hash(vode2)

    @test vode == similar(vode, t₀, q₀, p₀, λ₀)
    @test vode == similar(vode, t₀, q₀, p₀)
    @test vode == similar(vode, q₀, p₀)

end


@testset "$(rpad("Differential Algebraic Equations (DAE)",80))" begin

    dae_eqs = (dae_v, dae_u, dae_ϕ)

    dae  = DAE(dae_eqs..., t₀, [x₀], [λ₀])
    dae1 = DAE(dae_eqs..., t₀, x₀, λ₀)
    dae2 = DAE(dae_eqs..., x₀, λ₀)

    @test ndims(dae) == 2
    @test nsamples(dae) == 1
    @test nconstraints(dae) == 1
    @test periodicity(dae) == zero(x₀)
    @test get_function_tuple(dae) == NamedTuple{(:v, :u, :ϕ, :v̄)}((dae_eqs..., dae_v))

    @test dae == dae1
    @test dae == dae2

    @test hash(dae1) == hash(dae2)

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)

end


@testset "$(rpad("Partitioned Differential Algebraic Equations (PDAE)",80))" begin

    pdae_eqs = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    pdae  = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    pdae1 = PDAE(pdae_eqs..., t₀, q₀, p₀, λ₀)
    pdae2 = PDAE(pdae_eqs..., q₀, p₀, λ₀)

    @test ndims(pdae) == 1
    @test nsamples(pdae) == 1
    @test nconstraints(pdae) == 1
    @test periodicity(pdae) == zero(q₀)
    @test get_function_tuple(pdae) == NamedTuple{(:v, :f, :u, :g, :ϕ, :v̄, :f̄)}((pdae_eqs..., pdae_v, pdae_f))

    @test pdae == pdae1
    @test pdae == pdae2

    @test hash(pdae1) == hash(pdae2)

    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)

end


@testset "$(rpad("Implicit Differential Algebraic Equations (IDAE)",80))" begin

    idae_eqs = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    idae  = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v)
    idae1 = IDAE(idae_eqs..., t₀, q₀, p₀, λ₀; v̄=pdae_v)
    idae2 = IDAE(idae_eqs..., q₀, p₀, λ₀; v̄=pdae_v)

    @test ndims(idae) == 1
    @test nsamples(idae) == 1
    @test nconstraints(idae) == 1
    @test periodicity(idae) == zero(q₀)
    @test get_function_tuple(idae) == NamedTuple{(:ϑ, :f, :u, :g, :ϕ, :v̄, :f̄)}((idae_eqs..., pdae_v, pdae_f))

    @test idae == idae1
    @test idae == idae2

    @test hash(idae1) == hash(idae2)

    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)

end


@testset "$(rpad("Hamiltonian Differential Algebraic Equations (HDAE)",80))" begin

    hdae_eqs = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_u, pdae_g, pdae_ϕ, pdae_ψ, pdae_h)

    hdae  = HDAE(hdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    hdae1 = HDAE(hdae_eqs..., t₀, q₀, p₀, λ₀)
    hdae2 = HDAE(hdae_eqs..., q₀, p₀, λ₀)

    @test ndims(hdae) == 1
    @test nsamples(hdae) == 1
    @test nconstraints(hdae) == 1
    @test periodicity(hdae) == zero(q₀)
    @test get_function_tuple(hdae) == NamedTuple{(:v, :f, :u, :g, :u̅, :g̅, :ϕ, :ψ, :h, :v̄, :f̄)}((hdae_eqs..., pdae_v, pdae_f))

    @test hdae == hdae1
    @test hdae == hdae2

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)

end


@testset "$(rpad("Variational Differential Algebraic Equations (VDAE)",80))" begin

    vdae_eqs = (iode_ϑ, iode_f, iode_g, iode_g, pdae_ϕ, pdae_ψ)

    vdae  = VDAE(vdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
    vdae1 = VDAE(vdae_eqs..., t₀, q₀, p₀, λ₀, λ₀; v̄=iode_v)
    vdae2 = VDAE(vdae_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    vdae3 = VDAE(vdae_eqs..., t₀, q₀, p₀; v̄=iode_v)
    vdae4 = VDAE(vdae_eqs..., q₀, p₀; v̄=iode_v)

    @test ndims(vdae) == 1
    @test nsamples(vdae) == 1
    @test nconstraints(vdae) == 1
    @test periodicity(vdae) == zero(q₀)
    @test get_function_tuple(vdae) == NamedTuple{(:ϑ, :f, :g, :g̅, :ϕ, :ψ, :v̄, :f̄)}((vdae_eqs..., iode_v, iode_f))

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

end


@testset "$(rpad("Split Ordinary Differential Equations (SODE)",80))" begin

    f_sode = (f_sode_1, f_sode_2)
    q_sode = (q_sode_1, q_sode_2)

    sode  = SODE(f_sode, t₀, [q₀])
    sode1 = SODE(f_sode, t₀, q₀)
    sode2 = SODE(f_sode, q₀)

    @test ndims(sode) == 1
    @test nsamples(sode) == 1
    @test periodicity(sode) == zero(q₀)

    @test sode == sode1
    @test sode == sode2

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)

    @test sode == similar(sode, t₀, q₀)
    @test sode == similar(sode, q₀)


    sode  = SODE(f_sode, q_sode, t₀, [q₀])
    sode1 = SODE(f_sode, q_sode, t₀, q₀)
    sode2 = SODE(f_sode, q_sode, q₀)

    @test sode == sode1
    @test sode == sode2

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)

    @test sode == similar(sode, t₀, q₀)
    @test sode == similar(sode, q₀)
end


@testset "$(rpad("Split Partitioned Differential Algebraic Equations (SPDAE)",80))" begin

    spdae_eqs = ((pdae_v, pdae_u, pdae_u), (pdae_f, pdae_g, pdae_g), pdae_ϕ, pdae_ψ)

    spdae  = SPDAE(eltype(q₀), 1, 1, 1, 1, spdae_eqs..., t₀, q₀, p₀, λ₀)
    spdae1 = SPDAE(spdae_eqs..., t₀, q₀, p₀, λ₀)
    spdae2 = SPDAE(spdae_eqs..., t₀, q₀, p₀)
    spdae3 = SPDAE(spdae_eqs..., q₀, p₀)

    @test ndims(spdae) == 1
    @test nsamples(spdae) == 1
    @test nconstraints(spdae) == 1
    @test periodicity(spdae) == zero(q₀)
    @test get_function_tuple(spdae) == NamedTuple{(:v, :f, :ϕ, :ψ)}(spdae_eqs)

    @test spdae == spdae1
    @test spdae == spdae2
    @test spdae == spdae3

    @test hash(spdae) == hash(spdae1)
    @test hash(spdae) == hash(spdae2)
    @test hash(spdae) == hash(spdae3)

    @test spdae == similar(spdae, t₀, q₀, p₀, λ₀)
    @test spdae == similar(spdae, t₀, q₀, p₀)
    @test spdae == similar(spdae, q₀, p₀)

end

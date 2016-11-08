

t₀ = 0.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]


function f_ode(x, fx)
    fx[1] = x[1]
end

ode  = ODE{Float64}(1, 1, f_ode, t₀, q₀)
ode1 = ODE(f_ode, t₀, q₀)
ode2 = ODE(f_ode, q₀)

@test ode == ode1
@test ode == ode2


function f_pode(q, p, fp)
    fq[1] = q[1]
end

function g_pode(q, p, fp)
    fp[1] = 2p[1]
end

pode  = PODE{Float64}(1, 1, f_pode, g_pode, t₀, q₀, p₀)
pode1 = PODE(f_pode, g_pode, t₀, q₀, p₀)
pode2 = PODE(f_pode, g_pode, q₀, p₀)

@test pode == pode1
@test pode == pode2


sode  = SODE{Float64}(1, 1, q₀, p₀, t₀, do_nothing, do_nothing, f_pode, g_pode)
sode1 = SODE(q₀, p₀, t₀; v=f_pode, f=g_pode)
sode2 = SODE(q₀, p₀; v=f_pode, f=g_pode)

@test sode == sode1
@test sode == sode2
# @test pode1 == pode2


function f_dae(x, fx)
    fx[1] = x[1]
    fx[2] = x[2]
end

function u_dae(x, λ, fu)
    fu[1] = +λ[1]
    fu[1] = -λ[1]
end

function ϕ_dae(x, λ, fϕ)
    fϕ[1] = x[2] - x[1]
end

dae  = DAE{Float64}(2, 1, f_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
dae1 = DAE(f_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
dae2 = DAE(f_dae, u_dae, ϕ_dae, x₀, λ₀)

@test dae == dae1
@test dae == dae2


function f_pdae(q, p, fp)
    fq[1] = q[1]
end

function g_pdae(q, p, fp)
    fp[1] = p[1]
end

function u_pdae(q, p, λ, fu)
    fu[1] = +λ[1]
end

function v_pdae(q, p, λ, fv)
    fv[1] = -λ[1]
end

function ϕ_pdae(q, p, λ, fϕ)
    fϕ[1] = p[1] - q[1]
end

pdae  = PDAE{Float64}(1, 1, f_pdae, g_pdae, u_pdae, v_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
pdae1 = PDAE(f_pdae, g_pdae, u_pdae, v_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
pdae2 = PDAE(f_pdae, g_pdae, u_pdae, v_pdae, ϕ_pdae, q₀, p₀, λ₀)

@test pdae == pdae1
@test pdae == pdae2

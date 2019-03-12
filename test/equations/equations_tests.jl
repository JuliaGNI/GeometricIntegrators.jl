

t₀ = 0.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]


################################################################################
# Test ODE: Ordinary Differential Equation
################################################################################

function f_ode(t, x, f)
    f[1] = x[1]
end

ode  = ODE{eltype(q₀), typeof(t₀), typeof(f_ode), 1}(1, 1, f_ode, t₀, q₀)
ode1 = ODE(f_ode, t₀, q₀)
ode2 = ODE(f_ode, q₀)

@test ode == ode1
@test ode == ode2

@test hash(ode1) == hash(ode2)


################################################################################
# Test PODE: Partitioned Ordinary Differential Equation
################################################################################

function v_pode(t, q, p, v)
    v[1] = q[1]
end

function f_pode(t, q, p, f)
    f[1] = 2p[1]
end

pode  = PODE{eltype(q₀), typeof(t₀), typeof(v_pode), typeof(f_pode), 1}(1, 1, v_pode, f_pode, t₀, q₀, p₀)
pode1 = PODE(v_pode, f_pode, t₀, q₀, p₀)
pode2 = PODE(v_pode, f_pode, q₀, p₀)

@test pode == pode1
@test pode == pode2

@test hash(pode1) == hash(pode2)


################################################################################
# Test IODE: Implicit Ordinary Differential Equation
################################################################################

function iode_α(t, q, v, p)
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

iode  = IODE{eltype(q₀), typeof(t₀), typeof(iode_α), typeof(iode_f), typeof(iode_g), typeof(iode_v), 1}(1, 1, iode_α, iode_f, iode_g, iode_v, t₀, q₀, p₀, λ₀)
iode1 = IODE(iode_α, iode_f, iode_g, iode_v, t₀, q₀, p₀, λ₀)
iode2 = IODE(iode_α, iode_f, iode_g, iode_v, t₀, q₀, p₀)
iode3 = IODE(iode_α, iode_f, iode_g, iode_v, q₀, p₀)

@test iode == iode1
@test iode == iode2
@test iode == iode3

@test hash(iode1) == hash(iode2)


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

dae  = DAE{eltype(q₀), typeof(t₀), typeof(v_dae), typeof(u_dae), typeof(ϕ_dae), 1}(2, 1, 1, v_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
dae1 = DAE(v_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
dae2 = DAE(v_dae, u_dae, ϕ_dae, x₀, λ₀)

@test dae == dae1
@test dae == dae2

@test hash(dae1) == hash(dae2)


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
    u[1] = +λ[1]
end

function g_pdae(t, q, p, λ, g)
    g[1] = -λ[1]
end

function ϕ_pdae(t, q, p, λ, ϕ)
    ϕ[1] = p[1] - q[1]
end

function λ_pdae(t, q, p, v)
    nothing
end

pdae  = PDAE{eltype(q₀), typeof(t₀), typeof(v_pdae), typeof(f_pdae), typeof(u_pdae), typeof(g_pdae), typeof(ϕ_pdae), 1}(1, 1, 1, v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
pdae1 = PDAE(v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
pdae2 = PDAE(v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, q₀, p₀, λ₀)

@test pdae == pdae1
@test pdae == pdae2

@test hash(pdae1) == hash(pdae2)


idae  = IDAE{eltype(q₀), typeof(t₀), typeof(f_pdae), typeof(p_pdae), typeof(u_pdae), typeof(g_pdae), typeof(ϕ_pdae), typeof(λ_pdae), 1}(1, 1, 1, f_pdae, p_pdae, u_pdae, g_pdae, ϕ_pdae, λ_pdae, t₀, q₀, p₀, λ₀)
idae1 = IDAE(f_pdae, p_pdae, u_pdae, g_pdae, ϕ_pdae, λ_pdae, t₀, q₀, p₀, λ₀)
idae2 = IDAE(f_pdae, p_pdae, u_pdae, g_pdae, ϕ_pdae, λ_pdae, q₀, p₀, λ₀)

@test idae == idae1
@test idae == idae2

@test hash(idae1) == hash(idae2)


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



t₀ = 0.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]


function f_ode(t, x, f)
    f[1] = x[1]
end

ode  = ODE{eltype(q₀), typeof(t₀), typeof(f_ode)}(1, 1, f_ode, t₀, q₀)
ode1 = ODE(f_ode, t₀, q₀)
ode2 = ODE(f_ode, q₀)

@test ode == ode1
@test ode == ode2

@test hash(ode1) == hash(ode2)


function v_pode(t, q, p, v)
    v[1] = q[1]
end

function f_pode(t, q, p, f)
    f[1] = 2p[1]
end

pode  = PODE{eltype(q₀), typeof(t₀), typeof(v_pode), typeof(f_pode)}(1, 1, v_pode, f_pode, t₀, q₀, p₀)
pode1 = PODE(v_pode, f_pode, t₀, q₀, p₀)
pode2 = PODE(v_pode, f_pode, q₀, p₀)

@test pode == pode1
@test pode == pode2

@test hash(pode1) == hash(pode2)


iode  = IODE{eltype(q₀), typeof(t₀), typeof(v_pode), typeof(f_pode)}(1, 1, v_pode, f_pode, t₀, q₀, p₀)
iode1 = IODE(v_pode, f_pode, t₀, q₀, p₀)
iode2 = IODE(v_pode, f_pode, q₀, p₀)

@test iode == iode1
@test iode == iode2

@test hash(iode1) == hash(iode2)


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

dae  = DAE{eltype(q₀), typeof(t₀), typeof(v_dae), typeof(u_dae), typeof(ϕ_dae)}(2, 1, 1, v_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
dae1 = DAE(v_dae, u_dae, ϕ_dae, t₀, x₀, λ₀)
dae2 = DAE(v_dae, u_dae, ϕ_dae, x₀, λ₀)

@test dae == dae1
@test dae == dae2

@test hash(dae1) == hash(dae2)


function v_pdae(t, q, p, v)
    v[1] = q[1]
end

function f_pdae(t, q, p, f)
    f[1] = p[1]
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

pdae  = PDAE{eltype(q₀), typeof(t₀), typeof(v_pdae), typeof(f_pdae), typeof(u_pdae), typeof(g_pdae), typeof(ϕ_pdae)}(1, 1, 1, v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
pdae1 = PDAE(v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, t₀, q₀, p₀, λ₀)
pdae2 = PDAE(v_pdae, f_pdae, u_pdae, g_pdae, ϕ_pdae, q₀, p₀, λ₀)

@test pdae == pdae1
@test pdae == pdae2

@test hash(pdae1) == hash(pdae2)

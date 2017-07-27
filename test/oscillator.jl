
k = 0.5
ω = √k

q₀=[0.5, 0.0]
p₀=[0.0, 0.0]
λ₀=[0.0, 0.0]

A = sqrt(q₀[2]^2 / k + q₀[1]^2)
ϕ = asin(q₀[1] / A)


function oscillator_ode_f(t, x, f)
    f[1] = x[2]
    f[2] = -k*x[1]
    nothing
end

function oscillator_ode(x₀=q₀)
    ODE(oscillator_ode_f, x₀)
end


function oscillator_pode_v(t, q, p, v)
    v[1] = p[1]
    nothing
end

function oscillator_pode_f(t, q, p, f)
    f[1] = -k*q[1]
    nothing
end

function oscillator_pode(q₀=[q₀[1]], p₀=[p₀[1]])
    PODE(oscillator_pode_v, oscillator_pode_f, q₀, p₀)
end


function oscillator_sode_v_1(t, q, v)
    v[1] = q[2]
    v[2] = 0
    nothing
end

function oscillator_sode_v_2(t, q, v)
    v[1] = 0
    v[2] = -k*q[1]
    nothing
end

function oscillator_sode(q₀=q₀)
    SODE((oscillator_sode_v_1, oscillator_sode_v_2), q₀)
end


function oscillator_iode_α(t, q, v, p)
    p[1] = q[2]
    p[2] = 0
    nothing
end

function oscillator_iode_f(t, q, v, f)
    f[1] = -k*q[1]
    f[2] = v[1] - q[2]
    nothing
end

function oscillator_iode_g(t, q, λ, g)
    g[1] = 0
    g[2] = λ[1]
    nothing
end

function oscillator_iode_v(t, q, p, v)
    v[1] = q[2]
    v[2] = -k*q[1]
    nothing
end

function oscillator_iode(q₀=q₀, p₀=p₀)
    IODE(oscillator_iode_α, oscillator_iode_f,
         oscillator_iode_g, oscillator_iode_v,
         q₀, p₀)
end


function oscillator_idae_u(t, q, p, λ, u)
    u[1] = λ[1]
    u[2] = λ[2]
    nothing
end

function oscillator_idae_g(t, q, p, λ, g)
    g[1] = 0
    g[2] = λ[1]
    nothing
end

function oscillator_idae_ϕ(t, q, p, ϕ)
    ϕ[1] = p[1] - q[2]
    ϕ[2] = p[2]
    nothing
end

function oscillator_idae(q₀=q₀, p₀=p₀, λ₀=λ₀)
    IDAE(oscillator_iode_f, oscillator_iode_α,
         oscillator_idae_u, oscillator_idae_g,
         oscillator_idae_ϕ, oscillator_iode_v,
         q₀, p₀, λ₀)
end

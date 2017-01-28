
function lotka_volterra_2d_ode_f_params(t, x, f, a1, b1, a2, b2)
    f[1] = x[1] * (a1 + b1*x[2])
    f[2] = x[2] * (a2 + b2*x[1])
    nothing
end

function lotka_volterra_2d_ode(x₀=[1., 1.]; a1=-2., b1=+1., a2=+1., b2=-1.)
    function lotka_volterra_2d_ode_f(t, x, f)
       lotka_volterra_2d_ode_f_params(t, x, f, a1, b1, a2, b2)
    end
    ODE(lotka_volterra_2d_ode_f, x₀)
end


function lotka_volterra_2d_iode_p(t, q, v, p)
    p[1] = q[2] + log(q[2]) / q[1]
    p[2] = q[1]
    nothing
end

function lotka_volterra_2d_iode_f(t, q, v, f)
    f[1] = v[2] - log(q[2]) * v[1] / q[1]^2 - 1 + 1 / q[1]
    f[2] = v[1] * (1 + q[1] * q[2]) / (q[1] * q[2]) - 1 + 2 / q[2]
    nothing
end

function lotka_volterra_2d_iode(q₀=[1., 1.], p₀=[1. + log(1.), 1.])
    IODE(lotka_volterra_2d_iode_f, lotka_volterra_2d_iode_p, q₀, p₀)
end


function lotka_volterra_2d_idae_u(t, q, p, λ, u)
    u[1] = λ[1]
    u[2] = λ[2]
    nothing
end

function lotka_volterra_2d_idae_g(t, q, p, λ, g)
    g[1] = λ[2] - λ[1] * log(q[2]) / q[1]^2
    g[2] = λ[1] * (1 + 1 / (q[1] * q[2]) )
    nothing
end

function lotka_volterra_2d_idae_ϕ(t, q, p, ϕ)
    ϕ[1] = p[1] - q[2] - log(q[2]) / q[1]
    ϕ[2] = p[2] - q[1]
    nothing
end

function lotka_volterra_2d_idae_v(t, q, p, v)
    v[1] = q[1] * (q[2] - 2)
    v[2] = q[2] * (1 - q[1])
    nothing
end

function lotka_volterra_2d_idae(q₀=[1., 1.], p₀=[1. + log(1.), 1.], λ₀ = [0.0, 0.0])
    IDAE(lotka_volterra_2d_iode_f, lotka_volterra_2d_iode_p,
         lotka_volterra_2d_idae_u, lotka_volterra_2d_idae_g,
         lotka_volterra_2d_idae_ϕ, lotka_volterra_2d_idae_v, q₀, p₀, λ₀)
end


module LotkaVolterra2d

    using GeometricIntegrators.Equations

    export α1, α2, energy
    export lotka_volterra_2d_ode, lotka_volterra_2d_iode, lotka_volterra_2d_idae

    const A1=1.0
    const A2=1.0
    const B1=1.0
    const B2=2.0

    const X0=1.0
    const Y0=1.0


    function α1(t, q)
        q[2] + log(q[2]) / q[1]
    end

    function α2(t, q)
        q[1]
    end


    const Q0=[X0, Y0]
    const P0=[α1(0., Q0), α2(0., Q0)]


    function energy(t, q)
        q[1] + q[2] - log(q[1]) - 2*log(q[2])
    end


    function lotka_volterra_2d_ode_v(t, q, v)
        v[1] = q[1] * (A2*q[2] - B2)
        v[2] = q[2] * (B1 - A1*q[1])
        nothing
    end

    function lotka_volterra_2d_ode(q₀=Q0)
        ODE(lotka_volterra_2d_ode_v, q₀)
    end


    function lotka_volterra_2d_iode_α(t, q, v, p)
        p[1] = α1(t,q)
        p[2] = α2(t,q)
        nothing
    end

    function lotka_volterra_2d_iode_f(t, q, v, f)
        f[1] = v[2] - v[1] * log(q[2]) / q[1]^2 - A1 + B1 / q[1]
        f[2] = v[1] * (1 + 1 / (q[1] * q[2])) - A2 + B2 / q[2]
        nothing
    end

    function lotka_volterra_2d_iode_g(t, q, λ, g)
        g[1] = λ[2] - λ[1] * log(q[2]) / q[1]^2
        g[2] = λ[1] * (1 + 1 / (q[1] * q[2]))
        nothing
    end

    function lotka_volterra_2d_iode_v(t, q, p, v)
        v[1] = q[1] * (A2*q[2] - B2)
        v[2] = q[2] * (B1 - A1*q[1])
        nothing
    end

    function lotka_volterra_2d_iode(q₀=Q0, p₀=P0)
        IODE(lotka_volterra_2d_iode_α, lotka_volterra_2d_iode_f,
             lotka_volterra_2d_iode_g, lotka_volterra_2d_iode_v,
             q₀, p₀)
    end

    function lotka_volterra_2d_idae_u(t, q, p, λ, u)
        u[1] = λ[1]
        u[2] = λ[2]
        nothing
    end

    function lotka_volterra_2d_idae_g(t, q, p, λ, g)
        g[1] = λ[2] - λ[1] * log(q[2]) / q[1]^2
        g[2] = λ[1] * (1 + 1 / (q[1] * q[2]))
        nothing
    end

    function lotka_volterra_2d_idae_ϕ(t, q, p, ϕ)
        ϕ[1] = p[1] - α1(t,q)
        ϕ[2] = p[2] - α2(t,q)
        nothing
    end

    function lotka_volterra_2d_idae(q₀=Q0, p₀=P0, λ₀=zeros(Q0))
        IDAE(lotka_volterra_2d_iode_f, lotka_volterra_2d_iode_α,
             lotka_volterra_2d_idae_u, lotka_volterra_2d_idae_g,
             lotka_volterra_2d_idae_ϕ, lotka_volterra_2d_iode_v,
             q₀, p₀, λ₀)
    end

end

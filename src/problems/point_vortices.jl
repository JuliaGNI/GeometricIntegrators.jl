module PointVortices

    using GeometricIntegrators.Equations

    export point_vortices_ode, point_vortices_iode,
           hamiltonian, angular_momentum, α1, α2, α3, α4

    const γ₁ = +0.1
    const γ₂ = +0.1
    const q₀ = [1., +0.1, 1., -0.1]


    function S{T}(x::T,y::T)
        one(T) + x^2 + y^2
    end

    function S2{T}(x::T,y::T)
        one(T) + 2x^2 + 2y^2
    end

    function dSdx{T}(x::T,y::T)
        2x
    end

    function dSdy{T}(x::T,y::T)
        2y
    end


    function hamiltonian(t,q)
        γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)
    end

    function angular_momentum(t,q)
        # γ₁ * (q[1]^2 + q[2]^2) * S(q[1],q[2]) +
        # γ₂ * (q[3]^2 + q[4]^2) * S(q[3],q[4])
        q[1] * α2(t, q) - q[2] * α1(t, q) +
        q[3] * α4(t, q) - q[4] * α3(t, q)
    end


    function α1(t, q)
        - γ₁ * q[2] * S(q[1], q[2]) / 2
    end

    function α2(t, q)
        + γ₁ * q[1] * S(q[1], q[2]) / 2
    end

    function α3(t, q)
        - γ₂ * q[4] * S(q[3], q[4]) / 2
    end

    function α4(t, q)
        + γ₂ * q[3] * S(q[3], q[4]) / 2
    end


    const p₀ = [α1(0, q₀),
                α2(0, q₀),
                α3(0, q₀),
                α4(0, q₀)]


    function f1(t, q, v)
        γ₁ * ( dSdx(q[1],q[2]) * (q[1] * v[2] - q[2] * v[1]) + v[2] * S(q[1], q[2]) ) / 2
    end

    function f2(t, q, v)
        γ₁ * ( dSdy(q[1],q[2]) * (q[1] * v[2] - q[2] * v[1]) - v[1] * S(q[1], q[2]) ) / 2
    end

    function f3(t, q, v)
        γ₂ * ( dSdx(q[3],q[4]) * (q[3] * v[4] - q[4] * v[3]) + v[4] * S(q[3], q[4]) ) / 2
    end

    function f4(t, q, v)
        γ₂ * ( dSdy(q[3],q[4]) * (q[3] * v[4] - q[4] * v[3]) - v[3] * S(q[3], q[4]) ) / 2
    end


    function dHd1(t, q)
        + γ₁ * γ₂ * dSdx(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) +
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[1] - q[3]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dHd2(t, q)
        + γ₁ * γ₂ * dSdy(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) +
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[2] - q[4]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dHd3(t, q)
        + γ₁ * γ₂ * dSdx(q[3],q[4]) * S(q[1],q[2]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) -
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[1] - q[3]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dHd4(t, q)
        + γ₁ * γ₂ * dSdy(q[3],q[4]) * S(q[1],q[2]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) -
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[2] - q[4]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end


    function point_vortices_ode_v(t, q, v)
        denominator1 = 1 / (γ₁ * S2(q[1], q[2]))
        denominator2 = 1 / (γ₂ * S2(q[3], q[4]))

        v[1] = - dHd2(t,q) / denominator1
        v[2] = + dHd1(t,q) / denominator1
        v[3] = - dHd4(t,q) / denominator2
        v[4] = + dHd3(t,q) / denominator2

        nothing
    end

    function point_vortices_ode(q₀=q₀)
        ODE(point_vortices_ode_v, q₀)
    end


    function point_vortices_iode_α(t, q, v, p)
        p[1] = α1(t,q)
        p[2] = α2(t,q)
        p[3] = α3(t,q)
        p[4] = α4(t,q)
        nothing
    end

    function point_vortices_iode_f(t, q, v, f)
        f[1] = f1(t,q,v) - dHd1(t,q)
        f[2] = f2(t,q,v) - dHd2(t,q)
        f[3] = f3(t,q,v) - dHd3(t,q)
        f[4] = f4(t,q,v) - dHd4(t,q)
        nothing
    end

    function point_vortices_iode_g(t, q, λ, g)
        g[1] = f1(t,q,λ)
        g[2] = f2(t,q,λ)
        g[3] = f3(t,q,λ)
        g[4] = f4(t,q,λ)
        nothing
    end

    function point_vortices_iode_v(t, q, p, v)
        point_vortices_ode_v(t, q, v)
    end

    function point_vortices_iode(q₀=q₀, p₀=p₀)
        IODE(point_vortices_iode_α, point_vortices_iode_f,
             point_vortices_iode_g, point_vortices_iode_v,
             q₀, p₀)
    end

end

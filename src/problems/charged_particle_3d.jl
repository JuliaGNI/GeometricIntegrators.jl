
using GeometricIntegrators.Equations


function α1(t, q)
    q[4] + A1(t,q)
end

function α2(t, q)
    q[5] + A2(t,q)
end

function α3(t, q)
    q[6] + A3(t,q)
end


function α(t, q, p)
    p[1] = α1(t,q)
    p[2] = α2(t,q)
    p[3] = α3(t,q)
    p[4] = zero(eltype(q₀))
    p[5] = zero(eltype(q₀))
    p[6] = zero(eltype(q₀))
    nothing
end


function ϕ₀(x)
    E₀*sin(2π*x)
end

function ϕ(t,q)
   ϕ₀(q[3])
end

function E1(t,q)
   return zero(eltype(q))
end

function E2(t,q)
   return zero(eltype(q))
end

function E3(t,q)
   return - 2π*E₀*cos(2π*q[3])
end


function hamiltonian(t,q)
    0.5 * (q[4]^2 + q[5]^2 + q[6]^2) + ϕ(t, q)
end

function angular_momentum(t,q)
    q[1] * α2(t,q) - q[2] * α1(t,q)
end


function charged_particle_3d_iode_α(t, q, v, p)
    α(t, q, p)
end

function charged_particle_3d_iode_f(t, q, v, f)
    f[1] = dA1d1(t,q) * v[1] + dA2d1(t,q) * v[2] + dA3d1(t,q) * v[3] + E1(t,q)
    f[2] = dA1d2(t,q) * v[1] + dA2d2(t,q) * v[2] + dA3d2(t,q) * v[3] + E2(t,q)
    f[3] = dA1d3(t,q) * v[1] + dA2d3(t,q) * v[2] + dA3d3(t,q) * v[3] + E3(t,q)
    f[4] = v[1] - q[4]
    f[5] = v[2] - q[5]
    f[6] = v[3] - q[6]
    nothing
end

function charged_particle_3d_iode_g(t, q, λ, g)
    g[1] = dA1d1(t,q) * λ[1] + dA2d1(t,q) * λ[2] + dA3d1(t,q) * λ[3]
    g[2] = dA1d2(t,q) * λ[1] + dA2d2(t,q) * λ[2] + dA3d2(t,q) * λ[3]
    g[3] = dA1d3(t,q) * λ[1] + dA2d3(t,q) * λ[2] + dA3d3(t,q) * λ[3]
    g[4] = λ[1]
    g[5] = λ[2]
    g[6] = λ[3]
    nothing
end

function charged_particle_3d_iode_v(t, q, p, v)
    v[1] = q[4]
    v[2] = q[5]
    v[3] = q[6]
    v[4] = E1(t,q) + q[5] * B3(t,q) - q[6] * B2(t,q)
    v[5] = E2(t,q) + q[6] * B1(t,q) - q[4] * B3(t,q)
    v[6] = E3(t,q) + q[4] * B2(t,q) - q[5] * B1(t,q)
    nothing
end

function charged_particle_3d_iode(q₀=q₀)
    p₀ = zeros(q₀)

    if ndims(q₀) == 1
        α(0, q₀, p₀)
    else
        for i in 1:size(q₀,2)
            tq = zeros(eltype(q₀), size(q₀,1))
            tp = zeros(eltype(p₀), size(p₀,1))
            simd_copy_xy_first!(tq, q₀, i)
            α(0, tq, tp)
            simd_copy_yx_first!(tp, p₀, i)
        end
    end

    IODE(charged_particle_3d_iode_α, charged_particle_3d_iode_f,
         charged_particle_3d_iode_g, charged_particle_3d_iode_v,
         q₀, p₀)
end

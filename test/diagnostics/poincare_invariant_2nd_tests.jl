
module PoincareInvariant2ndTest

    using GeometricIntegrators
    using GeometricIntegrators.Utils
    using SymPy

    const nx = 100
    const ny = 100
    const Δt = 10.
    const B₀ = 1.
    const r₀ = 0.5
    const z₀ = 0.0
    const z₁ = 0.1
    const u₀ = 5E-1
    const u₁ = 5E-2
    const kz = 2
    const ku = 4


    function B(t, q)
        B₀ * (1 + q[1]^2 + q[2]^2)
    end

    function gcϑ(t, q, ϑ)
        ϑ[1] =-q[2] * (1 + B(t,q)) / 4
        ϑ[2] =+q[1] * (1 + B(t,q)) / 4
        ϑ[3] = q[4]
        ϑ[4] = 0
        nothing
    end

    function gcω(t, q, Ω)
        Ω[1,1] = 0
        Ω[1,2] =-B(t,q)
        Ω[1,3] = 0
        Ω[1,4] = 0

        Ω[2,1] =+B(t,q)
        Ω[2,2] = 0
        Ω[2,3] = 0
        Ω[2,4] = 0

        Ω[3,1] = 0
        Ω[3,2] = 0
        Ω[3,3] = 0
        Ω[3,4] =+q[4]

        Ω[4,1] = 0
        Ω[4,2] = 0
        Ω[4,3] =-q[4]
        Ω[4,4] = 0

        nothing
    end

    function gc_surface_q(s,t)
        x  = r₀*(s-0.5)
        y  = r₀*(t-0.5)
        z  = z₀ + z₁ * cos(2π*kz*s) * cos(2π*kz*t)
        u  = u₀ + u₁ * sin(2π*ku*s) * sin(2π*ku*t)

        [x, y, z, u]
    end

    function gc_surface_p(s,t)
        q = gc_surface_q(s,t)
        p = zeros(q)
        gcϑ(zero(eltype(q)), q, p)
        p
    end

    function gc_dummy(t, a, b, c)
        nothing
    end

    function gc_dummy_pode(q₀, p₀)
        PODE(gc_dummy, gc_dummy, q₀, p₀)
    end

    function gc_dummy_iode(q₀)
        p₀ = zeros(q₀)

        if ndims(q₀) == 1
            gcϑ(zero(eltype(q₀)), q₀, p₀)
        else
            tq = zeros(eltype(q₀), size(q₀,1))
            tp = zeros(eltype(p₀), size(p₀,1))

            for i in 1:size(q₀,2)
                simd_copy_xy_first!(tq, q₀, i)
                gcϑ(zero(eltype(q₀)), tq, tp)
                simd_copy_yx_first!(tp, p₀, i)
            end
        end

        IODE(gc_dummy, gc_dummy, gc_dummy, gc_dummy, q₀, p₀)
    end


    function compute_canonical_invariant_numerical()
        pinv = PoincareInvariant2ndCanonical(gc_dummy_pode, gc_surface_q, gc_surface_p, Δt, 4, nx, ny, 0)
        sol  = Solution(pinv.equ, Δt, 0)

        I = evaluate_poincare_invariant(pinv, sol)

        return I[1]
    end


    function compute_invariant_numerical()
        pinv = PoincareInvariant2nd(gc_dummy_iode, gc_surface_q, gcω, Δt, 4, nx, ny, 0)
        sol  = Solution(pinv.equ, Δt, 0)

        I, J, K, L = evaluate_poincare_invariant(pinv, sol)

        return I[1], J[1]
    end


    function compute_invariant_analytical()
        s, t = Sym("s, t")
        x, y, z, u = gc_surface_q(s,t)

        S = simplify(B₀ * (1 + x^2 + y^2) * ( diff(x, s) * diff(y, t) + diff(x, t) * diff(y, s) )
                   + u * ( diff(z, s) * diff(u, t) + diff(z, t) * diff(u, s) ) )
        I = N(SymPy.integrate(S, (s, 0, 1), (t, 0, 1)), 20)
    end


    export compute_canonical_invariant_numerical, compute_invariant_numerical, compute_invariant_analytical

end


using PoincareInvariant2ndTest


### compute analytical invariant ###

Ia = compute_invariant_analytical()


### compute and check canonical invariant ###

Jn = compute_canonical_invariant_numerical()

@assert Jn ≈ Ia atol=2eps()


### compute and check noncanonical invariant ###

In, Jn = compute_invariant_numerical()

@assert In ≈ Ia atol=2eps()
@assert Jn ≈ Ia atol=2eps()

# Which formula does each symplectic Euler method compute in one step? The expected values come from
# a Newton solve of the formulas, not from the package.
#
#   Euler-A:  q₁ = q₀ + h H_p(q₁, p₀),  p₁ = p₀ - h H_q(q₁, p₀)   Leimkuhler & Reich (2004), (4.8)-(4.9)
#   Euler-B:  q₁ = q₀ + h H_p(q₀, p₁),  p₁ = p₀ - h H_q(q₀, p₁)   Leimkuhler & Reich (2004), (4.6)-(4.7)
#
# Euler-B is the left and Euler-A the right form of Hairer, Lubich & Wanner (2006), Eq. (VI.3.1).
#
# Run with `julia --startup-file=no --project=. scripts/symplectic_euler_points.jl`. It exits
# nonzero when `IPRK` with a tableau misses its letter or reaches the other one, on a separable
# Hamiltonian when the explicit integrator does, and when `DVIA` or `DVIB` leaves the staggered form
# of its letter. On a non-separable Hamiltonian the explicit integrator computes neither formula:
# see `KNOWN_ISSUES.md`. The one-step map of `DVIA` and `DVIB` on (xₙ, yₙ) is no formula there
# either, and on a separable Hamiltonian it is the adjoint letter, which the script also checks.
using GeometricIntegrators
using Printf

function newton(g, x)
    for _ in 1:100
        δ = 1e-7 * (1 + abs(x))
        dx = g(x) / ((g(x + δ) - g(x - δ)) / (2δ))
        x -= dx
        abs(dx) < 1e-15 && return x
    end
    error("Newton did not converge")
end

function euler_a(Hq, Hp, q₀, p₀, h)
    q₁ = newton(q -> q - q₀ - h * Hp(q, p₀), q₀)
    (q₁, p₀ - h * Hq(q₁, p₀))
end

function euler_b(Hq, Hp, q₀, p₀, h)
    p₁ = newton(p -> p - p₀ + h * Hq(q₀, p), p₀)
    (q₀ + h * Hp(q₀, p₁), p₁)
end

function one_step(method, Hq, Hp, q₀, p₀, h)
    v!(v, t, q, p, params) = (v[1] = Hp(q[1], p[1]); nothing)
    f!(f, t, q, p, params) = (f[1] = -Hq(q[1], p[1]); nothing)
    sol = integrate(PODEProblem(v!, f!, (0.0, h), h, [q₀], [p₀]), method)
    (sol.q[end][1], sol.p[end][1])
end

dist(a, b) = hypot(a[1] - b[1], a[2] - b[2])

function check(label, x, A, B, letter)
    dA, dB = dist(x, A), dist(x, B)
    @printf("  %-28s %12.3e %12.3e\n", label, dA, dB)
    hit, other = letter === :A ? (dA, dB) : (dB, dA)
    (hit < 1e-12 && other > 1e-3) || error("$(label) is not Euler-$(letter)")
end

q₀, p₀, h = 0.7, 1.1, 0.1

# non-separable: H = ½p²(1 + q²) + cos q
Hq, Hp = (q, p) -> q * p^2 - sin(q), (q, p) -> p * (1 + q^2)
A, B = euler_a(Hq, Hp, q₀, p₀, h), euler_b(Hq, Hp, q₀, p₀, h)
println("non-separable H, distance to Euler-A and Euler-B")
check("IPRK(SymplecticEulerARK())",
    one_step(IPRK(SymplecticEulerARK()), Hq, Hp, q₀, p₀, h), A, B, :A)
check("IPRK(SymplecticEulerBRK())",
    one_step(IPRK(SymplecticEulerBRK()), Hq, Hp, q₀, p₀, h), A, B, :B)
for (label, method) in (("SymplecticEulerARK()", SymplecticEulerARK()),
    ("SymplecticEulerBRK()", SymplecticEulerBRK()))
    x = one_step(method, Hq, Hp, q₀, p₀, h)
    @printf("  %-28s %12.3e %12.3e\n", label, dist(x, A), dist(x, B))
end

# separable: H = ½p² + cos q
Hq, Hp = (q, p) -> -sin(q), (q, p) -> p
A, B = euler_a(Hq, Hp, q₀, p₀, h), euler_b(Hq, Hp, q₀, p₀, h)
println("separable H, distance to Euler-A and Euler-B")
check("SymplecticEulerARK()", one_step(SymplecticEulerARK(), Hq, Hp, q₀, p₀, h), A, B, :A)
check("SymplecticEulerBRK()", one_step(SymplecticEulerBRK(), Hq, Hp, q₀, p₀, h), A, B, :B)

# n steps of L = q₂ q̇₁ - H(q₁, q₂) with ϑ = (q₂, 0), as the states (xₖ, yₖ) = (q₁, q₂)
function dvi_states(method, H, Hx, Hy, x₀, y₀, h, n)
    ϑ!(θ, t, q, v, params) = (θ[1] = q[2]; θ[2] = 0; nothing)
    f!(f, t, q, v, params) = (f[1] = -Hx(q[1], q[2]); f[2] = v[1] - Hy(q[1], q[2]); nothing)
    g!(g, t, q, v, λ, params) = (g[1] = 0; g[2] = λ[1]; nothing)
    ω!(ω, t, q, v, params) = (fill!(ω, 0); ω[1, 2] = -1; ω[2, 1] = 1; nothing)
    l(t, q, v, params) = q[2] * v[1] - H(q[1], q[2])
    v̄!(v, t, q, p, params) = (v[1] = Hy(q[1], q[2]); v[2] = -Hx(q[1], q[2]); nothing)
    prob = LODEProblem(ϑ!, f!, g!, ω!, l, (0.0, n * h), h, [x₀, y₀], [y₀, 0.0]; v̄ = v̄!)
    sol = integrate(prob, method)
    [(sol.q[k][1], sol.q[k][2]) for k in 0:n]
end

# the largest distance from the states `z` to `euler` started at `z[1]`
function staggered_distance(euler, z, Hx, Hy, h)
    s, d = z[1], 0.0
    for zₖ in z[2:end]
        s = euler(Hx, Hy, s..., h)
        d = max(d, dist(s, zₖ))
    end
    d
end

function check_dvi(H, Hx, Hy, separable)
    A, B = euler_a(Hx, Hy, q₀, p₀, h), euler_b(Hx, Hy, q₀, p₀, h)
    zA = dvi_states(DVIA(), H, Hx, Hy, q₀, p₀, h, 10)
    zB = dvi_states(DVIB(), H, Hx, Hy, q₀, p₀, h, 10)
    if separable
        check("DVIA()", zA[2], A, B, :B)
        check("DVIB()", zB[2], A, B, :A)
    else
        @printf("  %-28s %12.3e %12.3e\n", "DVIA()", dist(zA[2], A), dist(zA[2], B))
        @printf("  %-28s %12.3e %12.3e\n", "DVIB()", dist(zB[2], A), dist(zB[2], B))
    end

    # DVIA on (xₙ₋₁, yₙ) is Euler-A and DVIB on (xₙ, yₙ₋₁) is Euler-B
    dA = staggered_distance(euler_a, [(zA[k][1], zA[k + 1][2]) for k in 1:10], Hx, Hy, h)
    dB = staggered_distance(euler_b, [(zB[k + 1][1], zB[k][2]) for k in 1:10], Hx, Hy, h)
    @printf("  %-28s %12.3e\n", "DVIA() from staggered A", dA)
    @printf("  %-28s %12.3e\n", "DVIB() from staggered B", dB)
    max(dA, dB) < 1e-13 || error("DVIA or DVIB is not the staggered form")
end

println("separable H, DVI, distance to Euler-A and Euler-B")
check_dvi((x, y) -> y^2 / 2 + cos(x), (x, y) -> -sin(x), (x, y) -> y, true)
println("non-separable H, DVI, distance to Euler-A and Euler-B")
check_dvi((x, y) -> y^2 * (1 + x^2) / 2 + cos(x), (x, y) -> x * y^2 - sin(x),
    (x, y) -> y * (1 + x^2), false)

println("COMPLETED")

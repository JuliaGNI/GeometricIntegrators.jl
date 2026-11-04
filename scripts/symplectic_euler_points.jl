# Which formula does each symplectic Euler method compute in one step? The expected values come from
# a Newton solve of the formulas, not from the package.
#
#   Euler-A:  q₁ = q₀ + h H_p(q₁, p₀),  p₁ = p₀ - h H_q(q₁, p₀)   Leimkuhler & Reich (2004), (4.8)-(4.9)
#   Euler-B:  q₁ = q₀ + h H_p(q₀, p₁),  p₁ = p₀ - h H_q(q₀, p₁)   Leimkuhler & Reich (2004), (4.6)-(4.7)
#
# Euler-B is the left and Euler-A the right form of Hairer, Lubich & Wanner (2006), Eq. (VI.3.1).
#
# Run with `julia --startup-file=no --project=scripts scripts/symplectic_euler_points.jl`. It exits
# nonzero when `IPRK` with a tableau misses its letter or reaches the other one, and on a separable
# Hamiltonian when the explicit integrator does. On a non-separable Hamiltonian the explicit
# integrator computes neither formula: see `KNOWN_ISSUES.md`.
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

println("COMPLETED")

# Symplecticity of the VSPARK projection methods — scope study.
#
#   julia --project=test scripts/vspark_projection_symplecticity.jl [steps...]
#
# Reproduces the S17 numbers of VERIFICATION_REPORT.md.
#
# The question is whether the projective-SPARK methods preserve the *noncanonical*
# symplectic structure of a degenerate Lagrangian system. Their theorem (SPARK Methods
# for Degenerate Lagrangian Systems, sec:projective_spark) requires four tableau
# conditions, b³ = b¹, b⁴ = b², AND ϕ(Q̃ᵢ,P̃ᵢ) = 0 at *every* projective stage. That
# manuscript states that none of its constructions satisfies all of them, so any
# symplecticity these methods have is not covered by the current theory.
#
# What is measured
# ----------------
# The first Poincaré invariant ∮p·dq on a closed loop of initial conditions, before and
# after n steps. For any symplectic map on (q,p) this is exactly conserved. Two variants:
#
#   canonical     ∮p·dq   with the integrator's own p    — the invariant of the map
#   noncanonical  ∮ϑ(q)·dq                               — the invariant of the
#                                                          constrained system
#
# The two agree whenever the method keeps p = ϑ(q), and separate when it does not.
#
# On the constraint manifold {p = ϑ(q)} the two coincide, so they agree to round-off for
# any method that preserves ϕ; they diverge for one that does not.
#
# The decisive diagnostic is the h dependence of the ONE-step defect: an exactly
# symplectic map has a defect at round-off independent of h, whereas a method with an
# O(hᵏ) defect shows clean decay. Step-count growth alone cannot distinguish the two.
#
# Two pitfalls this script avoids
# -------------------------------
# 1. dq/ds along the loop is differentiated *spectrally*. A central-difference stencil
#    leaves an O(N⁻²) quadrature error of order 1e-7 that masquerades as a defect.
# 2. Each problem is run in two gauges whose one-forms differ by an exact form. They
#    describe the same continuous system and share ∮ϑ·dq, but give different discrete
#    methods — and some methods are usable in one gauge and not the other.
#
# Optionally cross-checked against PoincareInvariants.jl if that package is available
# in the active environment (it is not a dependency of this repository).

using LinearAlgebra
using Printf

using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems
using GeometricProblems: LotkaVolterra2d, LotkaVolterra2dSingular,
    MasslessChargedParticle, MasslessChargedParticleSingular
using RungeKutta
import GeometricIntegratorsBase

const HAVE_PI = try
    @eval using PoincareInvariants
    true
catch
    false
end


# ---------------------------------------------------------------- problems ----
# Each entry supplies: an IDAEProblem built from a single q₀ (p₀ = ϑ(q₀) internally or
# explicitly), the one-form ϑ, the two-form Ω, and a loop of initial conditions.
#
# LotkaVolterra2d          ϑ = (q₂ + log q₂/q₁, q₁)
# LotkaVolterra2dSingular  ϑ = (log q₂/q₁, 0)                  differ by d(q₁q₂)
# MasslessChargedParticle          A = A₀/2 (1+x₁²+x₂²) (-x₂, x₁)
# MasslessChargedParticleSingular  A = -A₀ x₂ (1+2x₁²+⅔x₂²) (1, 0)   same B(x)

struct Problem
    name::String
    idae::Function      # q -> IDAEProblem with timespan/timestep
    theta::Function     # (t, q) -> ϑ
    omega::Function     # (t, q) -> Ω
    centre::Vector{Float64}
    radius::Float64
end

const LV_PARAMS = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
const MCP_PARAMS = (A₀ = 1.0, E₀ = 1.0)

function lv_problem(name, M)
    Problem(name,
        (q, T, h) -> M.idaeproblem(collect(q), M.ϑ(0.0, collect(q)), zero(collect(q));
            timespan = (0.0, T), timestep = h, parameters = LV_PARAMS),
        (t, q) -> M.ϑ(t, collect(q)),
        (t, q) -> (Ω = zeros(2, 2); M.lotka_volterra_2d_ω(Ω, t, collect(q), LV_PARAMS); Ω),
        [1.0, 1.0], 0.1)
end

function mcp_problem(name, M)
    Problem(name,
        (q, T, h) -> M.idaeproblem(collect(q); timespan = (0.0, T), timestep = h,
            parameters = MCP_PARAMS),
        (t, q) -> M.ϑ(t, collect(q), MCP_PARAMS),
        (t, q) -> M.ω(t, collect(q), MCP_PARAMS),
        [1.0, 1.0], 0.1)
end

const PROBLEMS = [
    lv_problem("LotkaVolterra2d", LotkaVolterra2d),
    lv_problem("LotkaVolterra2dSingular", LotkaVolterra2dSingular),
    mcp_problem("MasslessChargedParticle", MasslessChargedParticle),
    mcp_problem("MasslessChargedParticleSingular", MasslessChargedParticleSingular),
]


# ---------------------------------------------------------------- methods ----
# All projection variants that accept a Gauss inner method, at s = 1, 2, 3, plus the
# Lobatto-inner ones for contrast.

const GLRK_PROJECTIONS = [
    ("pSymplectic", TableauVSPARKGLRKpSymplectic),
    ("pSymmetric", TableauVSPARKGLRKpSymmetric),
    ("pMidpoint", TableauVSPARKGLRKpMidpoint),
    ("pModifiedMidpoint", TableauVSPARKGLRKpModifiedMidpoint),
    ("pInternal", TableauVSPARKGLRKpInternal),
    ("pModifiedInternal", TableauVSPARKGLRKpModifiedInternal),
    ("pLobattoIIIAIIIB", TableauVSPARKGLRKpLobattoIIIAIIIB),
    ("pLobattoIIIBIIIA", TableauVSPARKGLRKpLobattoIIIBIIIA),
    ("pModifiedLobattoIIIAIIIB", TableauVSPARKGLRKpModifiedLobattoIIIAIIIB),
    ("pModifiedLobattoIIIBIIIA", TableauVSPARKGLRKpModifiedLobattoIIIBIIIA),
]

const LOBATTO_PROJECTIONS = [
    ("LobIIIAIIIB pSymmetric", TableauVSPARKLobattoIIIAIIIBpSymmetric),
    ("LobIIIBIIIA pSymmetric", TableauVSPARKLobattoIIIBIIIApSymmetric),
    ("LobIIIAIIIB pMidpoint", TableauVSPARKLobattoIIIAIIIBpMidpoint),
    ("LobIIIBIIIA pMidpoint", TableauVSPARKLobattoIIIBIIIApMidpoint),
    ("LobIIIAIIIB pLobIIIAIIIB", TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB),
    ("LobIIIBIIIA pLobIIIAIIIB", TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB),
]

const N_LOOP = 128
const HS = (0.5, 0.2, 0.1, 0.05)

header(s) = println("\n", "="^110, "\n", s, "\n", "="^110)
sub(s) = println("\n--- ", s, " ", "-"^max(0, 104 - length(s)))


# ------------------------------------------------------- loop machinery ----

"Fourier spectral differentiation matrix for N uniform points on [0,2π)."
function fourier_diffmatrix(N)
    h = 2π / N
    D = zeros(N, N)
    for j in 1:N, k in 1:N
        j == k && continue
        D[j, k] = 0.5 * (-1)^(j + k) / tan((j - k) * h / 2)
    end
    D
end

"""
    loop_integral(qs, ps, D)

∮ p·dq on a closed loop sampled at `N` uniform values of the parameter s ∈ [0,1),
with dq/ds taken spectrally. For a smooth loop this is accurate to round-off.
"""
function loop_integral(qs, ps, D)
    N = length(qs)
    acc = 0.0
    for c in eachindex(qs[1])
        dqc = 2π .* (D * [q[c] for q in qs])
        acc += sum(ps[i][c] * dqc[i] for i in 1:N) / N
    end
    acc
end

circle(c, r, n) = [[c[1] + r * cos(2π * (i - 1) / n), c[2] + r * sin(2π * (i - 1) / n)]
                   for i in 1:n]

"""
    advance(prob, method, nsteps, h, D)

Integrate every point of the loop for `nsteps` steps and return the relative drift of
the canonical and noncanonical first invariants together with the largest constraint
violation `max|p - ϑ(q)|` on the image loop.
"""
function advance(prob::Problem, method, nsteps, h, D)
    qs = circle(prob.centre, prob.radius, N_LOOP)
    ps = [prob.theta(0.0, q) for q in qs]
    I0c = loop_integral(qs, ps, D)
    I0n = I0c                       # p = ϑ(q) on the initial loop

    q1 = Vector{Vector{Float64}}(undef, N_LOOP)
    p1 = Vector{Vector{Float64}}(undef, N_LOOP)
    ϕmax = 0.0
    for i in 1:N_LOOP
        sol = integrate(prob.idae(qs[i], nsteps * h, h), method)
        q1[i] = collect(sol.q[end])
        p1[i] = collect(sol.p[end])
        all(isfinite, q1[i]) && all(isfinite, p1[i]) || error("non-finite solution")
        ϕmax = max(ϕmax, maximum(abs, p1[i] .- prob.theta(nsteps * h, q1[i])))
    end
    I1c = loop_integral(q1, p1, D)
    I1n = loop_integral(q1, [prob.theta(nsteps * h, q) for q in q1], D)

    (can = abs(I1c - I0c) / abs(I0c), non = abs(I1n - I0n) / abs(I0n), ϕ = ϕmax)
end

"classify a one-step h sweep: round-off (flat) vs a genuine O(hᵏ) defect"
function classify(vals)
    ok = filter(isfinite, vals)
    isempty(ok) && return "fails"
    m = maximum(ok)
    m < 1e-13 && return length(ok) == length(vals) ? "round-off" : "round-off*"
    # fit the slope over the finite values that are above round-off
    length(ok) < 2 && return "defect"
    r = log2(maximum(ok) / max(minimum(ok), 1e-300)) / (length(ok) - 1) / log2(HS[1] / HS[2])
    @sprintf("defect O(h^%.0f)", max(r, 0.0))
end


# ---------------------------------------------------------------- step 1 ----
# The theorem's conditions, evaluated on the tableaus. Problem independent apart from
# the per-projective-stage requirement, which is measured on a converged step.

function tableau_conditions(m)
    t = GeometricIntegrators.SPARK.tableau(m)
    b1, b2 = t.p.b, t.p.β
    b3, b4 = t.q.b, t.q.β
    a1, a2 = t.p.a, t.p.α
    a3, a4 = t.q.a, t.q.α
    ã1, ã2 = t.p̃.a, t.p̃.α
    ã3, ã4 = t.q̃.a, t.q̃.α
    s, σ = length(b1), length(b2)
    (c1 = maximum(abs, [b1[i] * b3[j] - b1[i] * a3[i, j] - b3[j] * a1[j, i] for i in 1:s, j in 1:s]),
     c2 = maximum(abs, [b1[i] * b4[j] - b1[i] * a4[i, j] - b4[j] * ã1[j, i] for i in 1:s, j in 1:σ]),
     c3 = maximum(abs, [b2[i] * b3[j] - b2[i] * ã3[i, j] - b3[j] * a2[j, i] for i in 1:σ, j in 1:s]),
     c4 = maximum(abs, [b2[i] * b4[j] - b2[i] * ã4[i, j] - b4[j] * ã2[j, i] for i in 1:σ, j in 1:σ]),
     c5 = maximum(abs, b3 .- b1), c6 = maximum(abs, b4 .- b2))
end

function step1()
    header("Step 1 — the theorem's tableau conditions (problem independent)")
    println("  cond1..4 and b³=b¹, b⁴=b²; 0 means satisfied.")
    @printf("  %-34s %9s %9s %9s %9s %8s %8s\n",
        "method", "cond1", "cond2", "cond3", "cond4", "b³=b¹", "b⁴=b²")
    for (pname, ctor) in GLRK_PROJECTIONS, s in 1:3
        c = try
            tableau_conditions(ctor(s))
        catch
            nothing
        end
        c === nothing && continue
        @printf("  %-34s %9.1e %9.1e %9.1e %9.1e %8.0e %8.0e\n",
            "GLRK($s)$pname", c.c1, c.c2, c.c3, c.c4, c.c5, c.c6)
    end
    for (pname, ctor) in LOBATTO_PROJECTIONS, s in 2:3
        c = try
            tableau_conditions(ctor(s))
        catch
            nothing
        end
        c === nothing && continue
        @printf("  %-34s %9.1e %9.1e %9.1e %9.1e %8.0e %8.0e\n",
            "$pname($s)", c.c1, c.c2, c.c3, c.c4, c.c5, c.c6)
    end
    println()
    println("  Which condition fails depends on the projection, not only on R(∞):")
    println("    pSymplectic                cond4 fails for odd s (R(∞) = -1), holds at s = 2")
    println("    pSymmetric                 cond4 fails at every s (0.25) — its own docstring says so")
    println("    p{Modified,}Midpoint       cond2/cond3 fail for s = 1, 3; hold at s = 2")
    println("    pInternal                  cond2/cond3 fail at every s")
    println("    pModifiedInternal          cond2/cond3 fail for s = 1, 3; hold at s = 2")
    println("    p{Modified,}Lobatto*       all four hold at every s")
    println("  Only the per-projective-stage requirement ϕ(Q̃ᵢ,P̃ᵢ) = 0 then remains, and the")
    println("  ω of these constructions averages the projective constraints rather than")
    println("  imposing them individually.")
end


# ---------------------------------------------------------------- step 2 ----

function sweep(problems, methods, srange)
    D = fourier_diffmatrix(N_LOOP)
    for prob in problems
        sub(prob.name)
        @printf("  %-34s %9s %9s %9s %9s | %9s %9s %9s  %s\n",
            "method", "h=0.5", "h=0.2", "h=0.1", "h=0.05",
            "can 100", "non 100", "max|ϕ|", "verdict")
        for (pname, ctor) in methods, s in srange
            m = try
                ctor(s)
            catch
                continue
            end
            label = occursin(" ", pname) ? "$pname($s)" : "GLRK($s)$pname"
            vals = Float64[]
            @printf("  %-34s", label)
            for h in HS
                v = try
                    advance(prob, m, 1, h, D).can
                catch
                    NaN
                end
                push!(vals, v)
                @printf(" %9.2e", v)
            end
            r = try
                advance(prob, m, 100, 0.1, D)
            catch
                (can = NaN, non = NaN, ϕ = NaN)
            end
            @printf(" | %9.2e %9.2e %9.2e  %s\n", r.can, r.non, r.ϕ, classify(vals))
        end
    end
end

step2() = begin
    header("Step 2 — one-step canonical ∮p·dq defect vs h, Gauss inner method")
    println("  Flat at round-off across a factor of ten in h ⇒ exactly symplectic.")
    println("  Clean decay ⇒ a genuine O(hᵏ) defect. '100 st' is the 100-step drift at h=0.1.")
    sweep(PROBLEMS, GLRK_PROJECTIONS, 1:3)
end

step3() = begin
    header("Step 3 — the same, Lobatto inner methods (contrast)")
    sweep(PROBLEMS, LOBATTO_PROJECTIONS, 2:3)
end


# ---------------------------------------------------------------- step 4 ----
# Cross-check the inline spectral loop integral against PoincareInvariants.jl.

function step4()
    header("Step 4 — cross-check against PoincareInvariants.jl")
    if !HAVE_PI
        println("  PoincareInvariants.jl not available in this environment; skipped.")
        println("  To enable:  julia --project=test -e 'using Pkg; Pkg.add(\"PoincareInvariants\")'")
        return
    end
    D = fourier_diffmatrix(N_LOOP)
    prob = PROBLEMS[1]
    can = PoincareInvariants.FirstPI{Float64,4}(PoincareInvariants.canonical_one_form!, N_LOOP)
    qs = circle(prob.centre, prob.radius, N_LOOP)
    ps = [prob.theta(0.0, q) for q in qs]
    Z = Matrix{Float64}(undef, N_LOOP, 4)
    for i in 1:N_LOOP
        Z[i, 1:2] .= qs[i]
        Z[i, 3:4] .= ps[i]
    end
    mine = loop_integral(qs, ps, D)
    theirs = PoincareInvariants.compute!(can, Z, 0.0, nothing)
    @printf("  inline spectral   ∮p·dq = %.16e\n", mine)
    @printf("  PoincareInvariants ∮p·dq = %.16e\n", theirs)
    @printf("  relative difference       = %.2e\n", abs(mine - theirs) / abs(theirs))
end





# ---------------------------------------------------------------- step 5 ----
# Why do methods that violate the theorem's conditions still conserve the invariant?
#
# At the point where the manuscript invokes  P̃ᵢ = ϑ(Q̃ᵢ), the h¹ terms of the wedge
# product have been reduced to
#
#     Σᵢ b²ᵢ dG̃ᵢ ∧ dQ̃ᵢ + Σᵢ b⁴ᵢ dP̃ᵢ ∧ dΛ̃ᵢ
#   = Σᵢ (b⁴ᵢ - b²ᵢ) dϑ(Q̃ᵢ) ∧ dΛ̃ᵢ  +  Σᵢ b⁴ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ ,     Φ̃ᵢ = ϕ(Q̃ᵢ,P̃ᵢ) .
#
# With b⁴ = b² the first sum dies and what is actually required is the *two-form*
# condition
#
#     Σᵢ b⁴ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ = 0 .                                              (★)
#
# The theorem discharges (★) by the much stronger pointwise demand Φ̃ᵢ = 0. This step
# measures the two quantities that (★) is built from — the per-projective-stage
# constraint residual Φ̃ᵢ and the projective multiplier Λ̃ᵢ — so one can see which
# constructions satisfy (★) because a factor vanishes identically.

function projective_diagnostics(prob::Problem, m, h)
    q0 = prob.centre
    int = GeometricIntegrator(prob.idae(q0, h, h), m)
    ss = GeometricIntegratorsBase.solutionstep(int,
        GeometricIntegratorsBase.initialstate(prob.idae(q0, h, h)))
    GeometricIntegratorsBase.reset!(ss, h)
    GeometricIntegratorsBase.integrate!(ss, int)
    C = GeometricIntegratorsBase.cache(int)
    (Φ̃ = maximum(maximum(abs, x) for x in C.Φp),
     Λ̃ = maximum(maximum(abs, x) for x in C.Λp),
     Φ = maximum(maximum(abs, x) for x in C.Φi))
end

function step5()
    header("Step 5 — the two factors of the leftover term  Σᵢ b⁴ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ")
    println("  Φ̃ᵢ = ϕ(Q̃ᵢ,P̃ᵢ) at the projective stages (the theorem demands this be 0),")
    println("  Λ̃ᵢ = the projective multipliers. Either vanishing identically kills (★).")
    for prob in PROBLEMS
        sub(prob.name)
        @printf("  %-34s %11s %11s %11s\n", "method", "max|Φ̃ᵢ|", "max|Λ̃ᵢ|", "max|Φᵢ|")
        for (pname, ctor) in GLRK_PROJECTIONS, s in 1:3
            m = try ctor(s) catch; continue end
            r = try
                projective_diagnostics(prob, m, 0.1)
            catch
                (Φ̃ = NaN, Λ̃ = NaN, Φ = NaN)
            end
            @printf("  %-34s %11.2e %11.2e %11.2e\n", "GLRK($s)$pname", r.Φ̃, r.Λ̃, r.Φ)
        end
    end
end


const STEPS = isempty(ARGS) ? ["1", "2", "3", "4", "5"] : ARGS
"1" in STEPS && step1()
"2" in STEPS && step2()
"3" in STEPS && step3()
"4" in STEPS && step4()
"5" in STEPS && step5()

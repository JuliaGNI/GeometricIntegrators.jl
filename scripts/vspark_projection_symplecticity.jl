# Symplecticity of the VSPARK projection methods — scope study.
#
#   julia --project=scripts scripts/vspark_projection_symplecticity.jl [steps...]
#
# On a fresh checkout run  julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
#
# Reproduces the S17 numbers of docs/src/audit.md.
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
# 1. dq/ds along the loop is differentiated *spectrally* (see `loop_invariants.jl`). A
#    central-difference stencil leaves an O(N⁻²) quadrature error of order 1e-7 that
#    masquerades as a defect.
# 2. Each problem is run in two gauges whose one-forms differ by an exact form. They
#    describe the same continuous system and share ∮ϑ·dq, but give different discrete
#    methods — and some methods are usable in one gauge and not the other.
#
# All the numbers below come from the inline spectral loop integral. Step 4 additionally
# cross-checks that integral, on the *initial* loop only, against PoincareInvariants.jl —
# if and only if that package happens to be installed in the active environment. It is
# not a dependency of this repository and nothing here depends on it.

using LinearAlgebra
using Printf

using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems: LotkaVolterra2d, LotkaVolterra2dSingular,
    MasslessChargedParticle, MasslessChargedParticleSingular,
    PointVortices, PointVorticesLinear
using GeometricEquations: IDAEProblem
import GeometricIntegratorsBase

include(joinpath(@__DIR__, "loop_invariants.jl"))

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
        [1.0, 1.0], 0.1)
end

function mcp_problem(name, M)
    Problem(name,
        (q, T, h) -> M.idaeproblem(collect(q); timespan = (0.0, T), timestep = h,
            parameters = MCP_PARAMS),
        (t, q) -> M.ϑ(t, collect(q), MCP_PARAMS),
        [1.0, 1.0], 0.1)
end

# The point-vortex modules ship no IDAEProblem, so one is assembled here from their
# own ϑ, f and g plus the standard degenerate u and ϕ. PointVorticesLinear is the
# decisive case: its ϑ is LINEAR in q, which is the criterion below.
function pv_problem(name, M)
    prms = M.default_parameters()
    q0 = M.q₀
    idae = function (q, T, h)
        qq = collect(q)
        u = (u, t, q, v, p, λ, params) -> (u .= λ; nothing)
        ϕ = (ϕ, t, q, v, p, params) -> begin
            Θ = zero(p)
            M.point_vortices_ϑ(Θ, t, q, v, params)
            ϕ .= p .- Θ
            nothing
        end
        IDAEProblem(M.point_vortices_ϑ, M.point_vortices_f, u, M.point_vortices_g, ϕ,
            (0.0, T), h, qq, M.ϑ(qq), zero(qq);
            v̄ = M.point_vortices_v, parameters = prms)
    end
    Problem(name, idae, (t, q) -> M.ϑ(collect(q)), collect(q0), 0.1)
end

const PROBLEMS = [
    lv_problem("LotkaVolterra2d", LotkaVolterra2d),
    lv_problem("LotkaVolterra2dSingular", LotkaVolterra2dSingular),
    mcp_problem("MasslessChargedParticle", MasslessChargedParticle),
    mcp_problem("MasslessChargedParticleSingular", MasslessChargedParticleSingular),
    pv_problem("PointVorticesLinear  [ϑ LINEAR]", PointVorticesLinear),
    pv_problem("PointVortices", PointVortices),
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

"""
    selected_problems()

The entries of `PROBLEMS` picked out by a `--problems=a,b,...` argument; all of them by
default. Each pattern matches a problem name exactly if it can, otherwise as a
substring — so `--problems=MasslessChargedParticle` selects just that one and does not
also drag in `MasslessChargedParticleSingular`.
"""
function selected_problems()
    a = filter(startswith("--problems="), ARGS)
    isempty(a) && return PROBLEMS
    sel = Problem[]
    for pat in split(last(a)[12:end], ',')
        hit = filter(p -> p.name == pat, PROBLEMS)
        isempty(hit) && (hit = filter(p -> occursin(pat, p.name), PROBLEMS))
        isempty(hit) && error("--problems: `$pat` matches none of $(getproperty.(PROBLEMS, :name))")
        append!(sel, hit)
    end
    unique(sel)
end

header(s) = println("\n", "="^110, "\n", s, "\n", "="^110)
sub(s) = println("\n--- ", s, " ", "-"^max(0, 104 - length(s)))

# The singular and divergent cases below are the finding, not an accident, and their
# solver/line-search warnings would bury the tables. Each integration runs at
# `verbosity = 0` (merged into the method's `default_options`); the exception type is
# always printed instead, so nothing is hidden.


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
        sol = integrate(prob.idae(qs[i], nsteps * h, h), method; verbosity = 0, warn_iterations = 0)
        q1[i] = collect(sol.q[end])
        p1[i] = collect(sol.p[end])
        all(isfinite, q1[i]) && all(isfinite, p1[i]) || error("non-finite solution")
        ϕmax = max(ϕmax, maximum(abs, p1[i] .- prob.theta(nsteps * h, q1[i])))
    end
    I1c = loop_integral(q1, p1, D)
    I1n = loop_integral(q1, [prob.theta(nsteps * h, q) for q in q1], D)

    (can = abs(I1c - I0c) / abs(I0c), non = abs(I1n - I0n) / abs(I0n), ϕ = ϕmax)
end

# Anything at or below this is quadrature/solver round-off, not a defect. The loop
# integral is O(1e-2) here, so 1e-13 is ~1e-11 relative — a decade above eps.
const ROUNDOFF = 1e-13

"""
    classify(vals, hs = HS)

Classify a one-step `h` sweep: round-off (flat in `h`) versus a genuine `O(hᵏ)` defect.

`vals[i]` is the defect measured at `hs[i]`. Two effects have to be kept apart, and
getting either wrong invents an exponent:

  * `HS` is **not** geometrically uniform (0.5 → 0.2 is a factor 2.5, the rest are 2), so
    there is no single per-interval ratio to divide by; the exponent is fitted by least
    squares in `log h`.
  * A sweep can *bottom out*: the defect decays for the large `h` and then sits at
    round-off for the small ones. Including the floored points drags the fitted slope
    towards zero and reports, say, `O(h²)` for what is really a steep decay that ran out
    of dynamic range. Points at or below `ROUNDOFF` are therefore dropped from the fit
    and flagged with `†`, and a sweep that keeps fewer than two points is reported as
    `defect` with no exponent rather than a fabricated one.

A fitted slope that comes out non-positive means the defect does not decay with `h` at
all — the sweep is not measuring a convergence order, so it gets no exponent either,
rather than being clamped to a printable `O(h^0.0)`.

Suffixes: `*` some step size failed, `†` the fit used only the points above round-off.
"""
function classify(vals, hs = HS)
    @assert length(vals) == length(hs)
    finite = [i for i in eachindex(vals) if isfinite(vals[i])]
    isempty(finite) && return "fails"
    partial = length(finite) < length(vals) ? "*" : ""

    maximum(vals[i] for i in finite) <= ROUNDOFF && return "round-off" * partial

    keep = [i for i in finite if vals[i] > ROUNDOFF]
    floored = length(keep) < length(finite) ? "†" : ""
    length(keep) < 2 && return "defect" * floored * partial

    x = [log(hs[i]) for i in keep]
    y = [log(vals[i]) for i in keep]
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    r = sum((x .- x̄) .* (y .- ȳ)) / sum(abs2, x .- x̄)
    r > 0 || return "defect (no decay)" * floored * partial
    @sprintf("defect O(h^%.1f)%s%s", r, floored, partial)
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
            catch e
                @printf("  %-34s  no such tableau (%s)\n",
                    occursin(" ", pname) ? "$pname($s)" : "GLRK($s)$pname", typeof(e))
                continue
            end
            label = occursin(" ", pname) ? "$pname($s)" : "GLRK($s)$pname"
            vals = Float64[]
            why = String[]
            @printf("  %-34s", label)
            for h in HS
                v = try
                    advance(prob, m, 1, h, D).can
                catch e
                    push!(why, "h=$h:$(typeof(e))")
                    NaN
                end
                push!(vals, v)
                @printf(" %9.2e", v)
            end
            r = try
                advance(prob, m, 100, 0.1, D)
            catch e
                push!(why, "100st:$(typeof(e))")
                (can = NaN, non = NaN, ϕ = NaN)
            end
            @printf(" | %9.2e %9.2e %9.2e  %s\n", r.can, r.non, r.ϕ, classify(vals))
            isempty(why) || println(" "^38, "failures: ", join(why, ", "))
        end
    end
end

step2() = begin
    header("Step 2 — one-step canonical ∮p·dq defect vs h, Gauss inner method")
    println("  Flat at round-off across a factor of ten in h ⇒ exactly symplectic.")
    println("  Clean decay ⇒ a genuine O(hᵏ) defect. '100 st' is the 100-step drift at h=0.1.")
    sweep(selected_problems(), GLRK_PROJECTIONS, 1:3)
end

step3() = begin
    header("Step 3 — the same, Lobatto inner methods (contrast)")
    sweep(selected_problems(), LOBATTO_PROJECTIONS, 2:3)
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
    int = GeometricIntegrator(prob.idae(q0, h, h), m; verbosity = 0, warn_iterations = 0)
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
    for prob in selected_problems()
        sub(prob.name)
        @printf("  %-34s %11s %11s %11s\n", "method", "max|Φ̃ᵢ|", "max|Λ̃ᵢ|", "max|Φᵢ|")
        for (pname, ctor) in GLRK_PROJECTIONS, s in 1:3
            m = try ctor(s) catch; continue end
            note = ""
            r = try
                projective_diagnostics(prob, m, 0.1)
            catch e
                note = "  ($(typeof(e)))"
                (Φ̃ = NaN, Λ̃ = NaN, Φ = NaN)
            end
            @printf("  %-34s %11.2e %11.2e %11.2e%s\n", "GLRK($s)$pname", r.Φ̃, r.Λ̃, r.Φ, note)
        end
    end
end




# ---------------------------------------------------------------- step 6 ----
# Is (★) really the leftover?  Measure it directly and compare with the observed defect.
#
# Restricted to the constraint manifold, the initial data is q₀ ∈ R^D and every stage
# quantity is a function of q₀ (with p₀ = ϑ(q₀)). For D = 2 a two-form is a single
# number, the coefficient of dq₀¹ ∧ dq₀². Writing J[f] for the Jacobian ∂f/∂q₀,
#
#     (dα ∧ dβ)₁₂ = Σ_μ ( ∂α^μ/∂q₀¹ ∂β^μ/∂q₀² - ∂α^μ/∂q₀² ∂β^μ/∂q₀¹ ) ,
#
# so the predicted leftover of the proof is
#
#     star = h · Σᵢ b⁴ᵢ (dΦ̃ᵢ ∧ dΛ̃ᵢ)₁₂
#
# and the quantity it is supposed to account for is
#
#     defect = (dp₁ ∧ dq₁)₁₂ - (dp₀ ∧ dq₀)₁₂ .
#
# If star ≈ defect the diagnosis is right and (★) is the condition to attack. If not,
# the h² terms from the violated tableau conditions dominate instead.

"stage quantities of one step as a function of q₀, for finite differencing"
function stage_state(prob::Problem, m, q0, h)
    pr = prob.idae(q0, h, h)
    int = GeometricIntegrator(pr, m; verbosity = 0, warn_iterations = 0)
    ss = GeometricIntegratorsBase.solutionstep(int, GeometricIntegratorsBase.initialstate(pr))
    GeometricIntegratorsBase.reset!(ss, h)
    GeometricIntegratorsBase.integrate!(ss, int)
    C = GeometricIntegratorsBase.cache(int)
    cur = GeometricIntegratorsBase.current(ss)
    (Φ̃ = [copy(x) for x in C.Φp], Λ̃ = [copy(x) for x in C.Λp],
     q = collect(cur.q), p = collect(cur.p))
end

"(dα ∧ dβ)₁₂ from the two Jacobian columns of α and β"
wedge12(dα1, dα2, dβ1, dβ2) = sum(dα1 .* dβ2 .- dα2 .* dβ1)

function star_vs_defect(prob::Problem, m, h; ε = 1e-5)
    q0 = prob.centre
    pert(k, sgn) = (q = copy(q0); q[k] += sgn * ε; q)
    sp = [stage_state(prob, m, pert(k, +1), h) for k in 1:2]
    sm = [stage_state(prob, m, pert(k, -1), h) for k in 1:2]

    dΦ̃ = [[(sp[k].Φ̃[i] .- sm[k].Φ̃[i]) ./ (2ε) for k in 1:2] for i in eachindex(sp[1].Φ̃)]
    dΛ̃ = [[(sp[k].Λ̃[i] .- sm[k].Λ̃[i]) ./ (2ε) for k in 1:2] for i in eachindex(sp[1].Λ̃)]
    dq1 = [(sp[k].q .- sm[k].q) ./ (2ε) for k in 1:2]
    dp1 = [(sp[k].p .- sm[k].p) ./ (2ε) for k in 1:2]

    b4 = GeometricIntegrators.SPARK.tableau(m).q.β
    star = h * sum(b4[i] * wedge12(dΦ̃[i][1], dΦ̃[i][2], dΛ̃[i][1], dΛ̃[i][2])
                   for i in eachindex(dΦ̃))

    # (dp₀ ∧ dq₀)₁₂ with p₀ = ϑ(q₀), computed by the same stencil
    dϑ = [(prob.theta(0.0, pert(k, +1)) .- prob.theta(0.0, pert(k, -1))) ./ (2ε) for k in 1:2]
    e1 = [1.0, 0.0]
    e2 = [0.0, 1.0]
    before = wedge12(dϑ[1], dϑ[2], e1, e2)
    after = wedge12(dp1[1], dp1[2], dq1[1], dq1[2])
    (star = star, defect = after - before, scale = abs(before))
end

# The stage quantities come out of a nonlinear solve, so the central differences have a
# noise floor of roughly (solver tolerance)/ε on top of their O(ε²) truncation error.
# Both `star` and `defect` are therefore swept over ε: a value that is stable across two
# decades of ε is real, one that moves with ε is stencil noise and no ratio computed from
# it means anything.
const STAR_EPSILONS = (1e-4, 1e-5, 1e-6)

function step6()
    header("Step 6 — is (★) = h Σᵢ b⁴ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ the leftover?  (h = 0.2)")
    println("  'defect' is the measured one-step change of dp∧dq on the constraint manifold.")
    println("  If star ≈ defect, (★) is the term to attack; if not, the h² terms dominate.")
    println("  Swept over ε = ", STAR_EPSILONS, "; 'spread' is the relative variation of")
    println("  star/defect over the sweep. Trust the ratio only where the spread is small.")
    sel = [("pSymplectic", TableauVSPARKGLRKpSymplectic),
        ("pSymmetric", TableauVSPARKGLRKpSymmetric),
        ("pMidpoint", TableauVSPARKGLRKpMidpoint),
        ("pInternal", TableauVSPARKGLRKpInternal),
        ("pLobattoIIIAIIIB", TableauVSPARKGLRKpLobattoIIIAIIIB)]
    for prob in selected_problems()
        length(prob.centre) == 2 || continue   # the two-form bookkeeping assumes D = 2
        sub(prob.name)
        @printf("  %-30s %13s %13s %10s %8s\n",
            "method", "star", "defect", "star/defect", "spread")
        for (pname, ctor) in sel, s in 1:3
            m = try ctor(s) catch; continue end
            ratios = Float64[]
            rs = map(STAR_EPSILONS) do ε
                try
                    star_vs_defect(prob, m, 0.2; ε = ε)
                catch e
                    (star = NaN, defect = NaN, scale = NaN)
                end
            end
            for r in rs
                abs(r.defect) < 1e-14 && continue
                push!(ratios, r.star / r.defect)
            end
            r = rs[2]                       # report the middle ε, sweep the rest
            ratio = isempty(ratios) ? NaN : r.star / r.defect
            spread = length(ratios) < 2 ? NaN :
                     (maximum(ratios) - minimum(ratios)) / max(abs(ratio), eps())
            @printf("  %-30s %13.3e %13.3e %10.3f %8.1e\n",
                "GLRK($s)$pname", r.star, r.defect, ratio, spread)
        end
    end
end




# ---------------------------------------------------------------- step 7 ----
# The criterion that actually explains the round-off results.
#
# If every component of ϑ is at most LINEAR in q, an unprojected variational
# Runge-Kutta method already lands on the constraint manifold {p = ϑ(q)}: with
# ϑ(q) = Cq + d and b̄ = b,
#
#     ϕ(qₙ₊₁,pₙ₊₁) = h Σᵢ bᵢ [ (Cᵀ - C) Vₙ,ᵢ - ∇H(Qₙ,ᵢ) ] ,
#
# which the stage equations annihilate. The projection is then inert, Λ̃ = 0, and every
# Λ̃-carrying term in the wedge product drops out — so *any* projection built on such a
# method is symplectic, whichever of the theorem's conditions it violates.
#
# This step verifies the premise directly: an unprojected VPRK on the IODE form.

function step7()
    header("Step 7 — does an unprojected variational RK preserve ϕ by itself?")
    println("  Round-off ⇒ ϑ is (at most) linear in q, the projection is inert, and every")
    println("  projection method built on this inner method is symplectic for free.")
    cases = [
        ("PointVorticesLinear  [ϑ LINEAR]", () -> PointVorticesLinear.iodeproblem(),
            q -> PointVorticesLinear.ϑ(collect(q))),
        ("PointVortices", () -> PointVortices.iodeproblem(),
            q -> PointVortices.ϑ(collect(q))),
        ("LotkaVolterra2d", () -> LotkaVolterra2d.iodeproblem(),
            q -> LotkaVolterra2d.ϑ(0.0, collect(q))),
        ("MasslessChargedParticle", () -> MasslessChargedParticle.iodeproblem(),
            q -> MasslessChargedParticle.ϑ(0.0, collect(q), MCP_PARAMS)),
        ("MasslessChargedParticleSingular", () -> MasslessChargedParticleSingular.iodeproblem(),
            q -> MasslessChargedParticleSingular.ϑ(0.0, collect(q), MCP_PARAMS)),
    ]
    @printf("  %-34s %-16s %13s %13s\n", "problem", "method", "max|ϕ| 100st", "max|ϕ| 1000st")
    for (name, mk, th) in cases, (mn, m) in [("VPRKGauss(2)", VPRKGauss(2)), ("VPRKGauss(3)", VPRKGauss(3))]
        row = String[]
        for n in (100, 1000)
            v = try
                pr = similar(mk(); timespan = (0.0, n * 0.1), timestep = 0.1)
                sol = integrate(pr, m; verbosity = 0, warn_iterations = 0)
                maximum(maximum(abs, collect(sol.p[i]) .- th(sol.q[i])) for i in eachindex(sol.t))
            catch
                NaN
            end
            push!(row, @sprintf("%13.2e", v))
        end
        @printf("  %-34s %-16s %s %s\n", name, mn, row[1], row[2])
    end
    println("\n  Note: none of the Lotka-Volterra or charged-particle one-forms is linear, so the")
    println("  round-off results for those problems in steps 2-3 are NOT explained by this")
    println("  criterion — see step 5, where the multiplier Λ̃ is nonzero throughout.")
end


# `--problems=<substring>` restricts steps 2, 3, 5 and 6 to the matching problems, so a
# single row of the S17 tables can be reproduced without paying for all six. Every other
# argument selects a step.
const STEPS = let a = filter(!startswith("--"), ARGS)
    isempty(a) ? ["1", "2", "3", "4", "5", "6", "7"] : a
end

for k in ("1", "2", "3", "4", "5", "6", "7")
    k in STEPS && getfield(Main, Symbol("step", k))()
end

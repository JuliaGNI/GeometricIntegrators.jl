# Numerical confirmation of the SLRK verification findings.
#
#   julia --project=scripts scripts/slrk_verification.jl [steps...]
#
# On a fresh checkout run  julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
#
# Step 1 tableau symplecticity conditions + ω-matrix equivalence (problem independent)
# Step 2 S10: Jacobian singularity of the double-counted μ term at Δt = 1, before and
#        after the fix (the "before" residual is reconstructed, see
#        `residual_with_s10_bug!` — no source revert needed)
# Step 3 reference trajectories (compare across the S10 fix)
# Step 4 discrete symplecticity JᵀΩJ = Ω, constraint and energy drift
# Step 5 the same JᵀΩJ harness applied to VPRK (an invalid control — see below)
# Step 6 Poincaré invariant ∮p·dq (the valid test), plus the VPRK control
# Step 7 size of the leftover Λ-block term
#
# Steps 2-7 are run for every entry of `PROBLEMS`. `LotkaVolterra2d` and
# `LotkaVolterra2dSingular` describe the *same* continuous system: their one-forms
# differ by the exact form d(q₁q₂), so they share Ω = ∇ϑᵀ - ∇ϑ and the same
# Euler–Lagrange equations, but they give different discrete methods.
#
# Step 6 takes the step counts from the trailing arguments, so the long-time table of
# docs/src/audit.md is reproduced by
#
#   julia --project=scripts scripts/slrk_verification.jl 6 --steps=1,10,100,1000,3000

using ForwardDiff
using LinearAlgebra
using Printf

using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems: LotkaVolterra2d, LotkaVolterra2dSingular
using RungeKutta

import GeometricIntegratorsBase
using GeometricIntegratorsBase: solutionstep, current, history, nlsolution, reset!

include(joinpath(@__DIR__, "loop_invariants.jl"))

const PROBLEMS = [
    ("LotkaVolterra2d", LotkaVolterra2d),
    ("LotkaVolterra2dSingular", LotkaVolterra2dSingular),
]

const CONSTRUCTORS = [
    ("SLRKLobattoIIIAB", SLRKLobattoIIIAB),
    ("SLRKLobattoIIIBA", SLRKLobattoIIIBA),
    ("SLRKLobattoIIICC̄", SLRKLobattoIIICC̄),
    ("SLRKLobattoIIIC̄C", SLRKLobattoIIIC̄C),
    ("SLRKLobattoIIID", SLRKLobattoIIID),
    ("SLRKLobattoIIIE", SLRKLobattoIIIE),
]

const PARAMS = (a₁ = 1.0, a₂ = 1.0, b₁ = -1.0, b₂ = -2.0)
const Q0 = [1.0, 1.0]

header(s) = println("\n", "="^92, "\n", s, "\n", "="^92)
subheader(s) = println("\n--- ", s, " ", "-"^max(0, 86 - length(s)))

# The singular and divergent cases below are the finding, not an accident, and their
# solver/line-search warnings would bury the tables. Each integration therefore runs
# silently (`verbosity = 0, warn_iterations = 0`); the exception type is always printed
# instead, so nothing is hidden.

# problem-module accessors
ldae(M, q, p, λ; kwargs...) = M.ldaeproblem_slrk(q, p, λ; parameters = PARAMS, kwargs...)
lode(M, q, p; kwargs...) = M.lodeproblem(q, p; parameters = PARAMS, kwargs...)
theta(M, t, q) = M.ϑ(t, collect(q))
omega(M, t, q) = (Ω = zeros(2, 2); M.lotka_volterra_2d_ω(Ω, t, collect(q), PARAMS); Ω)
energy(M, t, q) = M.hamiltonian(t, collect(q), PARAMS)


# ---------------------------------------------------------------- step 1 ----
# The manuscript's symplecticity conditions collapse, for q̃ = q and p̃ = p, to
#   b̄ = b   and   diag(b̄) a + (diag(b) ā)ᵀ = b̄ bᵀ .

function check_pair(a, b, ā, b̄)
    (weights = maximum(abs, b̄ .- b),
        symplecticity = maximum(abs, Diagonal(b̄) * a .+ (Diagonal(b) * ā)' .- b̄ * b'))
end

function step1()
    header("Step 1a — tableau symplecticity conditions (problem independent)")
    worst = 0.0
    for (name, ctor) in CONSTRUCTORS, s in 2:4
        m = ctor(s)
        r1 = check_pair(m.q.a, m.q.b, m.p.a, m.p.b)
        r2 = check_pair(m.q̃.a, m.q̃.b, m.p̃.a, m.p̃.b)
        worst = max(worst, r1.weights, r1.symplecticity, r2.weights, r2.symplecticity)
        @printf("  %-18s s=%d   |b̄-b| = %.2e   |diag(b̄)a + (diag(b)ā)ᵀ - b̄bᵀ| = %.2e\n",
            name, s, max(r1.weights, r2.weights), max(r1.symplecticity, r2.symplecticity))
    end
    @printf("  → worst deviation over all constructors and s: %.2e\n", worst)

    header("Step 1b — ω matrix encodes {IIIA-averaged Ψ = 0, ϕ(qₙ₊₁,pₙ₊₁) = 0}")
    for s in 2:4
        ω = GeometricIntegrators.SPARK.lobatto_ω_matrix(s)
        A = TableauLobattoIIIA(s).a[2:s, 1:s]
        target = vcat(hcat(A, zeros(s - 1)), hcat(zeros(s)', 1))
        gap = maximum(abs, target * nullspace(ω)) + maximum(abs, ω * nullspace(target))
        @printf("  s=%d   size(ω) = %s   rank(ω) = %d   rank(target) = %d   cross-residual = %.2e\n",
            s, size(ω), rank(ω), rank(target), gap)
        s == 2 && println("        ω = ", round.(ω; digits = 12))
    end
end


# ---------------------------------------------------------------- step 2 ----

# The Jacobian is differentiated exactly (ForwardDiff, the same dual numbers the
# integrator's own solver uses). A central-difference stencil would put a noise floor
# of ~1e-9 under it, which turns the *exactly* singular Δt = 1 case into a finite
# κ ≈ 1e11 and understates the finding.
"""
    residual_with_s10_bug!(b, x, sol, params, int, method)

The *pre*-S10 residual: the null-vector multiplier μ added to the momentum-stage (`Z`) row
as well as to the primary-constraint (`Φ`) row.

The "before" column of the S10 table cannot be measured from the shipped source, because
the bug is fixed there. Rather than asking the reader to revert `integrators_slrk.jl` by
hand, it is reconstructed: `GeometricIntegratorsBase.residual!` gives the corrected
residual and this adds the removed term straight back, which reproduces the old code
exactly. The `Z` row of stage `i`, component `k`, is `b[4*(D*(i-1)+k-1)+2]`, and μ is the
trailing `D` unknowns of `x` — so `length(x) = D*(4s+1)`.
"""
function residual_with_s10_bug!(b, x, sol, params, int, method)
    GeometricIntegratorsBase.residual!(b, x, sol, params, int)
    S = method.s
    D = length(x) ÷ (4S + 1)
    for i in 1:S, k in 1:D
        b[4*(D*(i-1)+k-1)+2] += x[4*D*S+k] * method.d[i] / method.p.b[i]
    end
    b
end

function residual_jacobian(M, method, Δt; bug = false)
    prob = ldae(M, copy(Q0), theta(M, 0.0, Q0), zero(Q0); timespan = (0.0, 10Δt), timestep = Δt)
    int = GeometricIntegrator(prob, method)
    solstep = solutionstep(int, GeometricIntegratorsBase.initialstate(prob))
    reset!(solstep, Δt)

    sol = current(solstep)
    params = GeometricIntegratorsBase.parameters(solstep)
    GeometricIntegratorsBase.initial_guess!(sol, history(solstep), params, int)

    x = copy(nlsolution(int))
    r! = bug ?
         (b, y) -> residual_with_s10_bug!(b, y, sol, params, int, method) :
         (b, y) -> GeometricIntegratorsBase.residual!(b, y, sol, params, int)
    ForwardDiff.jacobian(r!, similar(x), x)
end

"""
    mu_visibility(J, D)

The smallest singular value of the last `D` columns of `J` — the null-vector multiplier
μ — after projecting out the span of all the other columns.

`cond(J)` alone says the stage system is ill conditioned but not *why*. This says it: it
is the size of the μ direction that survives in the reduced system, and with the S10 bug
it is proportional to `(1-Δt)` and vanishes at `Δt = 1`, because the momentum-stage
contribution then cancels the constraint-row one exactly.

The range of the other columns is taken from their SVD rather than a plain `qr`: those
columns are themselves rank deficient (that is the whole point of the null vector), and an
unpivoted QR gives no guarantee that its leading columns span the range.
"""
function mu_visibility(J, D)
    n = size(J, 2)
    rest = J[:, 1:(n-D)]
    mu = J[:, (n-D+1):n]
    F = svd(rest)
    keep = F.S .> max(size(rest)...) * eps() * maximum(F.S)
    U = F.U[:, keep]
    minimum(svdvals(mu .- U * (U' * mu)))
end

function step2()
    header("Step 2 — S10: conditioning of the stage Jacobian vs Δt")
    println("  κ ~ 1/(1-Δt) is the S10 signature; after the fix it stays flat.")
    println("  σ is `mu_visibility`: the μ direction the reduced system still sees.")
    println("  With the bug it is ∝ (1-Δt) and vanishes at Δt = 1.")
    println("  `before` re-adds the removed Z-row term (`residual_with_s10_bug!`);")
    println("  `after` is the shipped source. Both rows of the report's table.")
    D = length(Q0)
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS[1:2], s in 2:3
            m = ctor(s)
            for (label, bug) in (("before", true), ("after", false))
                @printf("  %-18s s=%d %-6s :", name, s, label)
                for Δt in (0.1, 0.5, 0.9, 0.99, 0.999, 1.0)
                    J = residual_jacobian(M, m, Δt; bug)
                    @printf("  Δt=%-5g κ=%.2e σ=%.2e", Δt, cond(J), mu_visibility(J, D))
                end
                println()
            end
        end
    end
end


# ---------------------------------------------------------------- step 3 ----

# A diverging solve does not always throw — `SLRKLobattoIIIBA(3)` on the singular gauge
# returns NaN instead — so the result is checked for finiteness rather than trusted.
function trajectory(M, method, Δt; T = 1.0)
    sol = integrate(ldae(M, copy(Q0), theta(M, 0.0, Q0), zero(Q0);
        timespan = (0.0, T), timestep = Δt), method; verbosity = 0, warn_iterations = 0)
    all(all(isfinite, sol.q[i]) && all(isfinite, sol.p[i]) for i in eachindex(sol.t)) ||
        error("non-finite solution")
    sol
end

function step3()
    header("Step 3 — reference trajectories at Δt = 0.1 (compare across the C1 fix)")
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            try
                sol = trajectory(M, ctor(s), 0.1)
                @printf("  %-18s s=%d   q(T) = [%.16e, %.16e]\n", name, s, sol.q[end][1], sol.q[end][2])
            catch e
                @printf("  %-18s s=%d   FAILED: %s\n", name, s, typeof(e))
            end
        end
    end
end


# ---------------------------------------------------------------- step 4 ----
# On the constraint manifold p = ϑ(q) one step is a map qₙ ↦ qₙ₊₁; the manuscript's
# theorem claims it preserves Ω. Valid for SLRK (which stays on the manifold), NOT
# for variational integrators (which do not) — see step 5.

function onestep_q(M, method, q0, Δt)
    prob = ldae(M, collect(q0), theta(M, 0.0, q0), zero(collect(q0));
        timespan = (0.0, Δt), timestep = Δt)
    q1 = collect(integrate(prob, method; verbosity = 0, warn_iterations = 0).q[end])
    all(isfinite, q1) || error("non-finite solution")
    q1
end

function jacobian_defect(M, method, q0, h, ε)
    Ω0 = omega(M, 0.0, q0)
    q1 = onestep_q(M, method, q0, h)
    J = zeros(2, 2)
    for j in 1:2
        qp = copy(q0); qp[j] += ε
        qm = copy(q0); qm[j] -= ε
        J[:, j] .= (onestep_q(M, method, qp, h) .- onestep_q(M, method, qm, h)) ./ (2ε)
    end
    maximum(abs, J' * omega(M, h, q1) * J .- Ω0)
end

function step4()
    header("Step 4a — discrete symplecticity ‖JᵀΩ(qₙ₊₁)J - Ω(qₙ)‖, ε sweep at Δt = 0.1")
    println("  (a real defect is ε-independent; solver/FD noise would scale like 1/ε)")
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            @printf("  %-18s s=%d :", name, s)
            for ε in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
                v = try
                    jacobian_defect(M, ctor(s), Q0, 0.1, ε)
                catch e
                    @printf("  %.0e→FAIL(%s)", ε, typeof(e))
                    NaN
                end
                isnan(v) || @printf("  %.0e→%.2e", ε, v)
            end
            println()
        end
    end

    header("Step 4b — Δt scaling of the symplecticity defect (ε = 1e-3)")
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            @printf("  %-18s s=%d :", name, s)
            prev = NaN
            for h in (0.2, 0.1, 0.05, 0.025)
                dev = try
                    jacobian_defect(M, ctor(s), Q0, h, 1e-3)
                catch e
                    @printf("  h=%-6g FAIL(%s)", h, typeof(e))
                    NaN
                end
                isnan(dev) && continue
                r = isnan(prev) ? NaN : log2(prev / dev)
                @printf("  h=%-6g %.2e (p≈%s)", h, dev, isnan(r) ? "--" : @sprintf("%.1f", r))
                prev = dev
            end
            println()
        end
    end

    header("Step 4c — constraint drift and energy drift over T = 100, Δt = 0.1")
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            try
                sol = trajectory(M, ctor(s), 0.1; T = 100.0)
                ϕmax = maximum(maximum(abs, sol.p[i] .- theta(M, sol.t[i], sol.q[i])) for i in eachindex(sol.t))
                H0 = energy(M, 0.0, sol.q[0])
                Hmax = maximum(abs(energy(M, sol.t[i], sol.q[i]) - H0) for i in eachindex(sol.t))
                @printf("  %-18s s=%d   max|ϕ| = %.3e   max|ΔH| = %.3e\n", name, s, ϕmax, Hmax)
            catch e
                @printf("  %-18s s=%d   FAILED: %s\n", name, s, typeof(e))
            end
        end
    end
end


# ---------------------------------------------------------------- step 5 ----
# CONTROL ATTEMPT THAT DOES NOT WORK. Variational integrators are symplectic on
# (q,p) space but do NOT preserve ϕ = p - ϑ(q), so restricting the map to the
# constraint manifold is invalid for them and the q-space test reports a spurious
# defect. Kept as a record of why step 6 is the right test.

function onestep_q_lode(M, method, q0, Δt)
    prob = lode(M, collect(q0), theta(M, 0.0, q0); timespan = (0.0, Δt), timestep = Δt)
    collect(integrate(prob, method; verbosity = 0, warn_iterations = 0).q[end])
end

function step5()
    header("Step 5 — INVALID CONTROL: JᵀΩJ for VPRK (variational integrators leave {ϕ=0})")
    controls = [("VPRKGauss(2)", VPRKGauss(2)),
        ("VPRKLobattoIIIAIIIĀ(2)", VPRKLobattoIIIAIIIĀ(2)),
        ("VPRKLobattoIIIAIIIĀ(3)", VPRKLobattoIIIAIIIĀ(3)),
        ("VPRKLobattoIIIBIIIB̄(2)", VPRKLobattoIIIBIIIB̄(2))]
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, m) in controls
            @printf("  %-24s:", name)
            for h in (0.2, 0.1, 0.05, 0.025)
                v = try
                    Ω0 = omega(M, 0.0, Q0)
                    q1 = onestep_q_lode(M, m, Q0, h)
                    J = zeros(2, 2)
                    for j in 1:2
                        qp = copy(Q0); qp[j] += 1e-3
                        qm = copy(Q0); qm[j] -= 1e-3
                        J[:, j] .= (onestep_q_lode(M, m, qp, h) .- onestep_q_lode(M, m, qm, h)) ./ 2e-3
                    end
                    maximum(abs, J' * omega(M, h, q1) * J .- Ω0)
                catch e
                    @printf("  h=%-6g FAIL(%s)", h, typeof(e))
                    NaN
                end
                isnan(v) || @printf("  h=%-6g %.2e", h, v)
            end
            println()
        end
    end
    println("\n  (nonzero here does NOT mean non-symplectic — the restriction is invalid.)")
end


# ---------------------------------------------------------------- step 6 ----
# Poincaré invariant. For ANY symplectic map on (q,p) the loop integral ∮p·dq is
# conserved — no constraint manifold, no finite-difference Jacobian. Valid for SLRK
# and VPRK alike, so the variational integrators are a genuine control.
#
# Two invariants are reported, because they are the same quantity only on {ϕ = 0}:
#
#   canonical     ∮p·dq    with the integrator's own p  — what a symplectic map conserves
#   noncanonical  ∮ϑ(q)·dq                              — the invariant of the
#                                                          constrained system
#
# For SLRK they agree to every digit (it holds ϕ to 1e-15), so either one measures its
# symplecticity. For VPRK they do NOT: a variational integrator conserves the canonical
# one exactly and leaves {ϕ = 0}, so its noncanonical value merely oscillates with
# max|ϕ| and is *not* a symplecticity control. Compare like with like — the decisive
# comparison is SLRK canonical against VPRK canonical.

function step6(step_counts = (1, 10, 100))
    header("Step 6 — Poincaré invariant ∮p·dq (exact for any symplectic map)")
    println("  canonical = ∮p·dq (the integrator's own p), noncanonical = ∮ϑ(q)·dq.")
    println("  For VPRK only the canonical column is a symplecticity control.")
    r = 0.1
    n = 200
    qs = circle(Q0, r, n)
    D = fourier_diffmatrix(n)

    for (pname, M) in PROBLEMS
        ps = [theta(M, 0.0, q) for q in qs]
        I0 = loop_integral(qs, ps, D)
        subheader("$pname   (loop: $n points, radius $r,  ∮p·dq = $(round(I0; sigdigits=12)))")

        "relative drift of both invariants after `nsteps` steps of size `h`"
        function advance(stepper, nsteps, h)
            q1 = Vector{Vector{Float64}}(undef, n)
            p1 = Vector{Vector{Float64}}(undef, n)
            ϕmax = 0.0
            for i in 1:n
                q1[i], p1[i] = stepper(qs[i], ps[i], nsteps, h)
                ϕmax = max(ϕmax, maximum(abs, p1[i] .- theta(M, nsteps * h, q1[i])))
            end
            Ican = loop_integral(q1, p1, D)
            Inon = loop_integral(q1, [theta(M, nsteps * h, q) for q in q1], D)
            (can = abs(Ican - I0) / abs(I0), non = abs(Inon - I0) / abs(I0), ϕ = ϕmax)
        end

        failed = (can = NaN, non = NaN, ϕ = NaN)

        slrk_step(m) = (q, p, nsteps, h) -> begin
            sol = integrate(ldae(M, copy(q), copy(p), zero(q);
                timespan = (0.0, nsteps * h), timestep = h), m; verbosity = 0, warn_iterations = 0)
            q1, p1 = collect(sol.q[end]), collect(sol.p[end])
            all(isfinite, q1) && all(isfinite, p1) || error("non-finite solution")
            (q1, p1)
        end
        lode_step(m) = (q, p, nsteps, h) -> begin
            sol = integrate(lode(M, copy(q), copy(p);
                timespan = (0.0, nsteps * h), timestep = h), m; verbosity = 0, warn_iterations = 0)
            q1, p1 = collect(sol.q[end]), collect(sol.p[end])
            all(isfinite, q1) && all(isfinite, p1) || error("non-finite solution")
            (q1, p1)
        end

        println("  CONTROL — symplectic variational integrators (canonical / noncanonical):")
        for (name, m) in [("VPRKGauss(2)", VPRKGauss(2)),
            ("VPRKGauss(3)", VPRKGauss(3)),
            ("VPRKLobattoIIIAIIIĀ(2)", VPRKLobattoIIIAIIIĀ(2)),
            ("VPRKLobattoIIIAIIIĀ(3)", VPRKLobattoIIIAIIIĀ(3))]
            @printf("    %-24s:", name)
            for nsteps in step_counts
                v = try
                    advance(lode_step(m), nsteps, 0.1)
                catch e
                    @printf("  %d: FAIL(%s)", nsteps, typeof(e))
                    failed
                end
                isnan(v.can) || @printf("  %4d: %.2e/%.2e", nsteps, v.can, v.non)
            end
            println()
        end

        println("  SLRK — accumulation over steps (canonical / noncanonical, max|ϕ|):")
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            @printf("    %-18s s=%d :", name, s)
            for nsteps in step_counts
                v = try
                    advance(slrk_step(ctor(s)), nsteps, 0.1)
                catch e
                    @printf("  %d: FAIL(%s)", nsteps, typeof(e))
                    failed
                end
                isnan(v.can) || @printf("  %4d: %.2e/%.2e(ϕ%.0e)", nsteps, v.can, v.non, v.ϕ)
            end
            println()
        end

        println("  SLRK — h dependence of the canonical defect after a single step:")
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            @printf("    %-18s s=%d :", name, s)
            prev = NaN
            for h in (0.2, 0.1, 0.05, 0.025)
                v = try
                    advance(slrk_step(ctor(s)), 1, h).can
                catch e
                    @printf("  h=%-6g FAIL(%s)", h, typeof(e))
                    NaN
                end
                isnan(v) && continue
                rr = isnan(prev) ? NaN : log2(prev / v)
                @printf("  h=%-6g %.2e(p≈%s)", h, v, isnan(rr) ? "--" : @sprintf("%.1f", rr))
                prev = v
            end
            println()
        end
    end
end


# ---------------------------------------------------------------- step 7 ----
# The proof leaves over  h · Σ_j β_j (d_j / b̄_j) · dμ ∧ dΛ_j  in the multiplier
# block. Σ_i d_i V_i = 0 is imposed; there is no such condition on Λ.

function step7()
    header("Step 7 — size of the leftover Λ-block term")
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            m = ctor(s)
            @printf("  %-18s s=%d :", name, s)
            for h in (0.2, 0.1, 0.05)
                prob = ldae(M, copy(Q0), theta(M, 0.0, Q0), zero(Q0); timespan = (0.0, h), timestep = h)
                int = GeometricIntegrator(prob, m; verbosity = 0, warn_iterations = 0)
                try
                    solstep = solutionstep(int, GeometricIntegratorsBase.initialstate(prob))
                    reset!(solstep, h)
                    GeometricIntegratorsBase.integrate!(solstep, int)
                    C = GeometricIntegratorsBase.cache(int)
                    dV = maximum(abs, sum(m.d[i] .* C.Vp[i] for i in 1:s))
                    dΛ = maximum(abs, sum(m.d[i] .* C.Λp[i] for i in 1:s))
                    μ = maximum(abs, C.μ)
                    @printf("  h=%-5g Σd·V=%.0e Σd·Λ=%.2e |μ|=%.1e h|μ||Σd·Λ|=%.2e",
                        h, dV, dΛ, μ, h * μ * dΛ)
                catch e
                    @printf("  h=%-5g FAIL(%s)", h, typeof(e))
                end
            end
            println()
        end
    end
    println("\n  (Σd·V ≈ 0 is the imposed constraint; Σd·Λ ≠ 0 is the uncontrolled leftover.")
    println("   h|μ||Σd·Λ| is the predicted size of the surviving wedge term.)")
end


# `--steps=1,10,100,...` sets the step counts of step 6; anything else selects steps.
const STEP_COUNTS = let a = filter(startswith("--steps="), ARGS)
    isempty(a) ? (1, 10, 100) : Tuple(parse.(Int, split(last(a)[9:end], ',')))
end
const STEPS = let a = filter(!startswith("--"), ARGS)
    isempty(a) ? ["1", "2", "3", "4", "5", "6", "7"] : a
end

"1" in STEPS && step1()
"2" in STEPS && step2()
"3" in STEPS && step3()
"4" in STEPS && step4()
"5" in STEPS && step5()
"6" in STEPS && step6(STEP_COUNTS)
"7" in STEPS && step7()

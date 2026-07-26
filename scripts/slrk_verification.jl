# Numerical confirmation of the SLRK verification findings.
#
#   julia --project=test scripts/slrk_verification.jl [steps...]
#
# Step 1 tableau symplecticity conditions + ω-matrix equivalence (problem independent)
# Step 2 C1: Jacobian singularity of the double-counted μ term at Δt = 1
# Step 3 reference trajectories (compare across the C1 fix)
# Step 4 discrete symplecticity JᵀΩJ = Ω, constraint and energy drift
# Step 5 the same JᵀΩJ harness applied to VPRK (an invalid control — see below)
# Step 6 Poincaré invariant ∮p·dq (the valid test), plus the VPRK control
# Step 7 size of the leftover Λ-block term
#
# Steps 2-7 are run for every entry of `PROBLEMS`. `LotkaVolterra2d` and
# `LotkaVolterra2dSingular` describe the *same* continuous system: their one-forms
# differ by the exact form d(q₁q₂), so they share Ω = ∇ϑᵀ - ∇ϑ and the same
# Euler–Lagrange equations, but they give different discrete methods.

using LinearAlgebra
using Printf

using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems
using GeometricProblems: LotkaVolterra2d, LotkaVolterra2dSingular
using RungeKutta

import GeometricIntegratorsBase
using GeometricIntegratorsBase: solutionstep, current, history, nlsolution, reset!

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

function residual_jacobian(M, method, Δt)
    prob = ldae(M, copy(Q0), theta(M, 0.0, Q0), zero(Q0); timespan = (0.0, 10Δt), timestep = Δt)
    int = GeometricIntegrator(prob, method)
    solstep = solutionstep(int, GeometricIntegratorsBase.initialstate(prob))
    reset!(solstep, Δt)

    sol = current(solstep)
    params = GeometricIntegratorsBase.parameters(solstep)
    GeometricIntegratorsBase.initial_guess!(sol, history(solstep), params, int)

    x = copy(nlsolution(int))
    n = length(x)
    J = zeros(n, n)
    ε = 1e-7
    rp = zeros(n)
    rm = zeros(n)
    for j in 1:n
        xp = copy(x); xp[j] += ε
        xm = copy(x); xm[j] -= ε
        GeometricIntegratorsBase.residual!(rp, xp, sol, params, int)
        GeometricIntegratorsBase.residual!(rm, xm, sol, params, int)
        J[:, j] .= (rp .- rm) ./ (2ε)
    end
    J
end

function step2()
    header("Step 2 — C1: conditioning of the stage Jacobian vs Δt")
    for (pname, M) in PROBLEMS
        subheader(pname)
        for (name, ctor) in CONSTRUCTORS[1:2], s in 2:3
            @printf("  %-18s s=%d :", name, s)
            for Δt in (0.1, 0.5, 0.9, 0.99, 0.999, 1.0)
                @printf("  Δt=%-5g κ=%.2e", Δt, cond(residual_jacobian(M, ctor(s), Δt)))
            end
            println()
        end
    end
    println("\n  (κ ~ 1/(1-Δt) is the C1 signature; after the fix it stays flat.)")
end


# ---------------------------------------------------------------- step 3 ----

function trajectory(M, method, Δt; T = 1.0)
    integrate(ldae(M, copy(Q0), theta(M, 0.0, Q0), zero(Q0); timespan = (0.0, T), timestep = Δt), method)
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
    collect(integrate(prob, method).q[end])
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
                catch
                    NaN
                end
                @printf("  %.0e→%.2e", ε, v)
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
                catch
                    NaN
                end
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
    collect(integrate(prob, method).q[end])
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
                catch
                    NaN
                end
                @printf("  h=%-6g %.2e", h, v)
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
∮ p·dq on a closed loop sampled at N uniform values of the parameter s ∈ [0,1).

dq/ds is taken spectrally, so for a smooth loop the quadrature is accurate to
round-off; a central-difference stencil leaves an O(N⁻²) error of order 1e-7 that
masquerades as a symplecticity defect.
"""
function loop_integral(qs, ps, D = fourier_diffmatrix(length(qs)))
    N = length(qs)
    acc = 0.0
    for c in eachindex(qs[1])
        dqc = 2π .* (D * [q[c] for q in qs])
        acc += sum(ps[i][c] * dqc[i] for i in 1:N) / N
    end
    acc
end

circle(q0, r, n) = [[q0[1] + r * cos(2π * (i - 1) / n), q0[2] + r * sin(2π * (i - 1) / n)] for i in 1:n]

function step6()
    header("Step 6 — Poincaré invariant ∮p·dq (exact for any symplectic map)")
    r = 0.1
    n = 200
    qs = circle(Q0, r, n)
    D = fourier_diffmatrix(n)

    for (pname, M) in PROBLEMS
        ps = [theta(M, 0.0, q) for q in qs]
        I0 = loop_integral(qs, ps, D)
        subheader("$pname   (loop: $n points, radius $r,  ∮p·dq = $(round(I0; sigdigits=12)))")

        function advance(stepper, nsteps, h)
            q1 = Vector{Vector{Float64}}(undef, n)
            p1 = Vector{Vector{Float64}}(undef, n)
            for i in 1:n
                q1[i], p1[i] = stepper(qs[i], ps[i], nsteps, h)
            end
            loop_integral(q1, p1, D)
        end

        slrk_step(m) = (q, p, nsteps, h) -> begin
            sol = integrate(ldae(M, copy(q), copy(p), zero(q);
                    timespan = (0.0, nsteps * h), timestep = h), m)
            (collect(sol.q[end]), collect(sol.p[end]))
        end
        lode_step(m) = (q, p, nsteps, h) -> begin
            sol = integrate(lode(M, copy(q), copy(p);
                    timespan = (0.0, nsteps * h), timestep = h), m)
            (collect(sol.q[end]), collect(sol.p[end]))
        end

        println("  CONTROL — symplectic variational integrators:")
        for (name, m) in [("VPRKGauss(2)", VPRKGauss(2)),
            ("VPRKLobattoIIIAIIIĀ(2)", VPRKLobattoIIIAIIIĀ(2)),
            ("VPRKLobattoIIIAIIIĀ(3)", VPRKLobattoIIIAIIIĀ(3))]
            @printf("    %-24s:", name)
            for (nsteps, h) in ((1, 0.1), (10, 0.1), (100, 0.1))
                I1 = try
                    advance(lode_step(m), nsteps, h)
                catch
                    NaN
                end
                @printf("  %3d×%g: %.2e", nsteps, h, abs(I1 - I0) / abs(I0))
            end
            println()
        end

        println("  SLRK — accumulation over steps:")
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            @printf("    %-18s s=%d :", name, s)
            for (nsteps, h) in ((1, 0.1), (10, 0.1), (100, 0.1))
                I1 = try
                    advance(slrk_step(ctor(s)), nsteps, h)
                catch
                    NaN
                end
                @printf("  %3d×%g: %.2e", nsteps, h, abs(I1 - I0) / abs(I0))
            end
            println()
        end

        println("  SLRK — h dependence after a single step:")
        for (name, ctor) in CONSTRUCTORS, s in 2:3
            @printf("    %-18s s=%d :", name, s)
            prev = NaN
            for h in (0.2, 0.1, 0.05, 0.025)
                v = try
                    abs(advance(slrk_step(ctor(s)), 1, h) - I0) / abs(I0)
                catch
                    NaN
                end
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
                int = GeometricIntegrator(prob, m)
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


const STEPS = isempty(ARGS) ? ["1", "2", "3", "4", "5", "6", "7"] : ARGS

"1" in STEPS && step1()
"2" in STEPS && step2()
"3" in STEPS && step3()
"4" in STEPS && step4()
"5" in STEPS && step5()
"6" in STEPS && step6()
"7" in STEPS && step7()

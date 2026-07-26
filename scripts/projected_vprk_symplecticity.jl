#!/usr/bin/env julia
#
# Symplecticity of projected variational partitioned Runge-Kutta methods
# ======================================================================
#
# Numerical verification of the symplecticity properties of the projection
# methods for degenerate Lagrangian systems described in
#
#   M. Kraus, "Projected Variational Integrators for Degenerate Lagrangian
#   Systems", Sections "Standard Projection" ... "Internal Projection".
#
# For a degenerate Lagrangian  L(q,q̇) = ϑ(q)⋅q̇ - H(q)  the Dirac constraint is
#
#   φ(q,p) = p - ϑ(q) = 0 ,
#
# and the projected integrator is a map  Φ_h : Δ → Δ,  q_n ↦ q_{n+1}, where the
# momenta and the Lagrange multiplier are internal variables.  Symplecticity of
# Φ_h means preservation of the noncanonical two-form on M,
#
#   ω̄ = ½ Ω̄_ij(q) dq^i ∧ dq^j ,      Ω̄_ij = ∂ϑ_j/∂q^i - ∂ϑ_i/∂q^j ,
#
# i.e.  J^T Ω̄(q_{n+1}) J = Ω̄(q_n)  with  J = ∂q_{n+1}/∂q_n.  The exception is
# the symplectic projection, which is symplectic on the extended space M × R^d
# with respect to the modified form ω_λ (see `Omega_lambda` below).
#
# This script implements one step of each projected method directly from the
# equations of the paper -- independently of the integrators in
# GeometricIntegrators -- and evaluates the defect in the symplecticity
# condition together with its order in h.  The step map is differentiated
# exactly (implicit function theorem + ForwardDiff), so the reported defects are
# not polluted by finite differences; the harness is validated by checking that
# the *unprojected* VPRK method is canonically symplectic to round-off.
#
# Results obtained with this script (massless charged particle, GLRK1; the two
# Lotka-Volterra models are degenerate for this test, see `report_degeneracy`):
#
#   scheme       form            defect (h=0.2)  order   paper's claim    verdict
#   standard     ω̄                     7.2e-08      4    not symplectic   confirmed
#   symmetric    ω̄                     2.2e-10    ~4.5   pseudo-sympl.    confirmed
#   symplectic   dp∧dq - dA∧dB         1.9e-16     --    symplectic       confirmed (exact)
#   symplectic   ω_λ, mixed sign +     1.9e-16     --    --               exact after sign fix
#   symplectic   ω_λ as printed        6.7e-05      3    --               printed sign is wrong
#   midpoint     ω̄                     8.9e-06      3    SYMPLECTIC       REFUTED
#   midpoint_R1  --                 unsolvable     --    --               ill-posed
#   internal     ω̄                     8.9e-06      3    (empty section)  not symplectic
#
# The midpoint projection is the interesting case: its symplecticity proof needs
# R(∞) = +1, while the requirement that the midpoint coincide with an internal
# stage forces R(∞) = -1 (see `check_midpoint_stage_lemma` below).  The two
# hypotheses are mutually exclusive, and with the sign that makes the method
# solvable the symplecticity defect is O(h^3) per step.  Dropping the R(∞) factor
# instead (as printed in the arXiv version) makes the step equations unsolvable.
#
# For methods with a midpoint stage the multiplier λ -- and hence the defect -- is
# governed by the same midpoint/trapezoidal discrepancy irrespective of the order
# of the base method, which is why GLRK1 and SRK3 give nearly identical defects.
#
# Usage:  julia --project=scripts scripts/projected_vprk_symplecticity.jl

module ProjectedVPRKSymplecticity

using ForwardDiff
using LinearAlgebra
using Printf

using GeometricProblems.LotkaVolterra2d
using GeometricProblems.LotkaVolterra2dSingular
using GeometricProblems.MasslessChargedParticle


# ---------------------------------------------------------------------------
# degenerate Lagrangian test problems
# ---------------------------------------------------------------------------

struct DegenerateProblem{TT,TH}
    name::String
    ϑ::TT
    H::TH
    q₀::Vector{Float64}
end

function testproblems()
    lv, lvs, mcp = LotkaVolterra2d, LotkaVolterra2dSingular, MasslessChargedParticle
    par, pars, parm = lv.default_parameters(), lvs.default_parameters(), mcp.default_parameters()
    [
        DegenerateProblem("LotkaVolterra2d",
            q -> lv.ϑ(0.0, q), q -> lv.hamiltonian(0.0, q, par), copy(lv.q₀)),
        DegenerateProblem("LotkaVolterra2dSingular",
            q -> lvs.ϑ(0.0, q), q -> lvs.hamiltonian(0.0, q, pars), copy(lvs.q₀)),
        # Both Lotka-Volterra models turn out to be degenerate for this test (see
        # `report_degeneracy`), so the massless charged particle of the same paper
        # is included as a discriminating example: there both components of ϑ are
        # nonlinear.
        DegenerateProblem("MasslessChargedParticle",
            q -> mcp.ϑ(0.0, q, parm), q -> mcp.hamiltonian(0.0, q, parm), copy(mcp.q₀)),
    ]
end

# φ(q,p) = p - ϑ(q);  the projection direction is
# Ω^{-1} ∇φ^T λ = (φ_p^T λ, -φ_q^T λ) = (λ, ∇ϑ(q)^T λ).
# `fϑ(prob, q, v)` returns (∇ϑ(q)^T v)_i = ∂ϑ_k/∂q^i v^k, computed as a
# gradient so that it nests cheaply inside further differentiation.
fϑ(prob, q, v) = ForwardDiff.gradient(x -> dot(prob.ϑ(x), v), q)
∇H(prob, q) = ForwardDiff.gradient(prob.H, q)

# Ω̄_ij = ∂ϑ_j/∂q^i - ∂ϑ_i/∂q^j
function Omega_bar(prob, q)
    J = ForwardDiff.jacobian(prob.ϑ, q)
    transpose(J) - J
end

# Modified two-form of the symplectic projection on M × R^d.  The paper gives
#   ω_λ = ½ Ω̄_ij dq^i∧dq^j - (h²/2) Ω̄_ij dλ^i∧dλ^j - h² λ^k ϑ_{k,ij} dq^i∧dλ^j ,
# but the exactly conserved form dp∧dq - dA∧dB (see `check_generalised_form`) has
# the opposite sign on the mixed term.  With `mixedsign = +1` the form below is
# preserved to round-off; with the printed sign `mixedsign = -1` a spurious defect
# of order h³ appears, since the two differ by 2 h² λ^k ϑ_{k,ij} = O(h³).
function Omega_lambda(prob, q, λ, h; mixedsign=+1)
    Ω = Omega_bar(prob, q)
    T = ForwardDiff.hessian(x -> dot(prob.ϑ(x), λ), q)   # (λ⋅ϑ_qq)_ij = λ^k ϑ_{k,ij}
    [Ω (mixedsign * h^2 .* T); (-mixedsign * h^2 .* T) (-h^2 .* Ω)]
end


# ---------------------------------------------------------------------------
# tableaus of variational partitioned Runge-Kutta methods
# ---------------------------------------------------------------------------

struct VPRKTableau
    name::String
    a::Matrix{Float64}
    ā::Matrix{Float64}
    b::Vector{Float64}
    R∞::Float64
    midpointstage::Int      # index m with Z_m = ½(z̄_n + z̄_{n+1}), or 0
end

"""
Build a VPRK tableau from `(a, b)`.  The partner coefficients follow from the
symplecticity conditions  b_i ā_ij + b_j a_ji = b_i b_j,  b̄ = b, and
R(∞) = 1 - bᵀA⁻¹e is computed rather than tabulated.
"""
function VPRKTableau(name, a, b)
    s = length(b)
    ā = [b[j] * (1 - a[j, i] / b[i]) for i in 1:s, j in 1:s]
    R∞ = 1 - dot(b, a \ ones(s))
    m = 0
    for i in 1:s
        isapprox(a[i, :], b ./ 2; atol=1e-13) && (m = i)
    end
    VPRKTableau(name, a, ā, b, R∞, m)
end

function tableaus()
    s3 = sqrt(3) / 6
    s15 = sqrt(15)
    [
        VPRKTableau("GLRK1", [0.5;;], [1.0]),
        VPRKTableau("GLRK2", [1/4 1/4-s3; 1/4+s3 1/4], [0.5, 0.5]),
        VPRKTableau("GLRK3",
            [5/36 2/9-s15/15 5/36-s15/30
             5/36+s15/24 2/9 5/36-s15/24
             5/36+s15/30 2/9+s15/15 5/36], [5/18, 4/9, 5/18]),
        VPRKTableau("SRK3",                        # Table `tab:srk3` of the paper
            [5/36 2/9 25/180-s15/10
             5/36 2/9 5/36
             25/180+s15/10 2/9 5/36], [5/18, 4/9, 5/18]),
    ]
end


# ---------------------------------------------------------------------------
# projection schemes
# ---------------------------------------------------------------------------
#
# All endpoint projections share the structure
#
#   0  = φ(z_n)
#   z̄_n     = z_n     + h RU₁ Ω^{-1}∇φ^T(ζ₁) λ₁
#   z̄_{n+1} = Ψ_h(z̄_n)
#   z_{n+1} = z̄_{n+1} + h RU₂ Ω^{-1}∇φ^T(ζ₂) λ₂
#   0  = φ(z_{n+1})
#
# and differ only in the weights (RU₁, RU₂), the evaluation points (ζ₁, ζ₂) and
# whether the pre- and post-multipliers λ₁, λ₂ are the same variable.
#
#   :standard    RU = (0,  1 ),  ζ₂ = z_{n+1}
#   :symmetric   RU = (1, R∞),  ζ₁ = z_n, ζ₂ = z_{n+1},          λ₁ = λ₂ = λ_{n+1/2}
#   :symplectic  RU = (1, R∞),  ζ₁ = z_n, ζ₂ = z_{n+1},          λ₁ = λ_n, λ₂ = λ_{n+1}
#   :midpoint    RU = (1, R∞),  ζ₁ = ζ₂ = ½(z̄_n + z̄_{n+1}),      λ₁ = λ₂ = λ_{n+1/2}
#   :midpoint_R1 as :midpoint but with RU₂ = +1 (the sign the proof of the paper
#                requires, and the one printed in the arXiv version) -- unsolvable
#   :midpoint_gi as :midpoint but ζ = ½(z_n + z_{n+1}), i.e. the midpoint of the
#                *projected* points, as implemented in GeometricIntegrators
#
# The internal projection (`:internal`) perturbs the internal stages instead of
# the endpoints and is handled separately.

const ENDPOINT_SCHEMES = (:standard, :symmetric, :symplectic, :midpoint, :midpoint_R1, :midpoint_gi)
const ALL_SCHEMES = (ENDPOINT_SCHEMES..., :internal)

extended(scheme) = scheme === :symplectic          # map on M × R^d instead of M
claimed(scheme) = scheme === :symplectic || scheme === :midpoint || scheme === :midpoint_R1

function weights(scheme, R∞)
    scheme === :standard && return (0.0, 1.0)
    scheme === :midpoint_R1 && return (1.0, 1.0)
    return (1.0, R∞)
end

nunknowns(scheme, s, d) = scheme === :internal ? d * (s + 1) : d * (s + 2)

"""
    residual(x, qₙ, λₙ, h, tab, prob, scheme)

Residual of one step of the projected VPRK method.  Unknowns are
`x = [V₁…V_s; p̄ₙ; λ]` for the endpoint projections and `x = [V₁…V_s; λ]` for the
internal projection.  Returns `(residual, qₙ₊₁, λ)`.
"""
function residual(x, qₙ, λₙ, h, tab::VPRKTableau, prob, scheme::Symbol)
    d = length(qₙ)
    s = length(tab.b)
    V = [x[d*(i-1)+1:d*i] for i in 1:s]
    pₙ = prob.ϑ(qₙ)

    if scheme === :internal
        λ = x[d*s+1:d*s+d]
        # every internal stage is shifted by h Ω^{-1}∇φ^T λ
        Q = [qₙ .+ h .* sum(tab.a[i, j] .* V[j] for j in 1:s) .+ h .* λ for i in 1:s]
        F = [fϑ(prob, Q[i], V[i]) .- ∇H(prob, Q[i]) for i in 1:s]
        G = sum(tab.b[j] .* fϑ(prob, Q[j], λ) for j in 1:s)
        qₙ₊₁ = qₙ .+ h .* sum(tab.b[i] .* V[i] for i in 1:s) .+ (h * (1 + tab.R∞)) .* λ
        pₙ₊₁ = pₙ .+ h .* sum(tab.b[i] .* F[i] for i in 1:s) .+ (h * (1 + tab.R∞)) .* G
        res = vcat([prob.ϑ(Q[i]) .- pₙ .- h .* sum(tab.ā[i, j] .* F[j] for j in 1:s) .- h .* G for i in 1:s]...,
            pₙ₊₁ .- prob.ϑ(qₙ₊₁))
        return res, qₙ₊₁, λ
    end

    p̄ₙ = x[d*s+1:d*s+d]
    λ = x[d*s+d+1:d*s+2d]
    RU₁, RU₂ = weights(scheme, tab.R∞)
    λ₁ = scheme === :symplectic ? λₙ : λ            # pre-perturbation multiplier

    # perturbation (q-component needs no evaluation point since φ_p = 1)
    q̄ₙ = qₙ .+ (h * RU₁) .* λ₁

    # one step of the variational partitioned Runge-Kutta method
    Q = [q̄ₙ .+ h .* sum(tab.a[i, j] .* V[j] for j in 1:s) for i in 1:s]
    F = [fϑ(prob, Q[i], V[i]) .- ∇H(prob, Q[i]) for i in 1:s]
    q̄ₙ₊₁ = q̄ₙ .+ h .* sum(tab.b[i] .* V[i] for i in 1:s)
    p̄ₙ₊₁ = p̄ₙ .+ h .* sum(tab.b[i] .* F[i] for i in 1:s)

    # projection
    qₙ₊₁ = q̄ₙ₊₁ .+ (h * RU₂) .* λ
    ζ = if scheme === :midpoint || scheme === :midpoint_R1
        (q̄ₙ .+ q̄ₙ₊₁) ./ 2                          # midpoint of the perturbed points
    elseif scheme === :midpoint_gi
        (qₙ .+ qₙ₊₁) ./ 2                          # midpoint of the projected points
    else
        qₙ₊₁
    end
    ζ₁ = (scheme === :midpoint || scheme === :midpoint_R1 || scheme === :midpoint_gi) ? ζ : qₙ
    pₙ₊₁ = p̄ₙ₊₁ .+ (h * RU₂) .* fϑ(prob, ζ, λ)

    res = vcat([prob.ϑ(Q[i]) .- p̄ₙ .- h .* sum(tab.ā[i, j] .* F[j] for j in 1:s) for i in 1:s]...,
        p̄ₙ .- pₙ .- (h * RU₁) .* fϑ(prob, ζ₁, λ₁),   # perturbation, p-component
        pₙ₊₁ .- prob.ϑ(qₙ₊₁))                        # 0 = φ(z_{n+1})

    return res, qₙ₊₁, λ
end


# ---------------------------------------------------------------------------
# solving one step and differentiating the step map
# ---------------------------------------------------------------------------

"velocity of the underlying Euler-Lagrange equations, Ω̄(q) q̇ = ∇H(q)"
qdot(prob, q) = Omega_bar(prob, q) \ ∇H(prob, q)

function initialguess(prob, qₙ, h, tab, scheme)
    d = length(qₙ)
    s = length(tab.b)
    x = zeros(nunknowns(scheme, s, d))
    v = qdot(prob, qₙ)
    for i in 1:s
        x[d*(i-1)+1:d*i] .= v
    end
    scheme === :internal || (x[d*s+1:d*s+d] .= prob.ϑ(qₙ))
    x
end

"Newton solver with exact (ForwardDiff) Jacobian; returns `(x, ‖res‖, converged)`."
function newton(prob, qₙ, λₙ, h, tab, scheme; itmax=100, tol=1e-13)
    x = initialguess(prob, qₙ, h, tab, scheme)
    f = y -> residual(y, qₙ, λₙ, h, tab, prob, scheme)[1]
    nr = norm(f(x))
    for _ in 1:itmax
        nr < tol && break
        J = ForwardDiff.jacobian(f, x)
        Δ = try
            J \ f(x)
        catch
            return x, nr, false
        end
        all(isfinite, Δ) || return x, nr, false
        # damped step: the R(∞)=+1 midpoint variant has a singular Jacobian
        α = 1.0
        for _ in 1:20
            xt = x .- α .* Δ
            nt = norm(f(xt))
            if isfinite(nt) && nt < nr
                x, nr = xt, nt
                break
            end
            α /= 2
        end
        α < 1e-6 && break
    end
    x, nr, nr < 1e-10
end

"""
    stepmap(prob, qₙ, λₙ, h, tab, scheme)

Solve one step and return `(J, qₙ₊₁, λₙ₊₁, ‖res‖, converged)` where `J` is the
Jacobian of the step map: `∂qₙ₊₁/∂qₙ` for the projections that act on Δ, and
`∂(qₙ₊₁,λₙ₊₁)/∂(qₙ,λₙ)` for the symplectic projection.  The Jacobian is exact:
the implicit function theorem is applied to the converged solution.
"""
function stepmap(prob, qₙ, λₙ, h, tab, scheme)
    d = length(qₙ)
    x, nr, ok = newton(prob, qₙ, λₙ, h, tab, scheme)
    _, qout, λout = residual(x, qₙ, λₙ, h, tab, prob, scheme)
    ok || return nothing, qout, λout, nr, false

    out = extended(scheme) ? (r -> vcat(r[2], r[3])) : (r -> r[2])
    y = extended(scheme) ? vcat(qₙ, λₙ) : qₙ
    splity = extended(scheme) ? (v -> (v[1:d], v[d+1:2d])) : (v -> (v, λₙ))

    Fx = ForwardDiff.jacobian(z -> residual(z, qₙ, λₙ, h, tab, prob, scheme)[1], x)
    gx = ForwardDiff.jacobian(z -> out(residual(z, qₙ, λₙ, h, tab, prob, scheme)), x)
    Fy = ForwardDiff.jacobian(v -> begin
        q, l = splity(v)
        residual(x, q, l, h, tab, prob, scheme)[1]
    end, y)
    gy = ForwardDiff.jacobian(v -> begin
        q, l = splity(v)
        out(residual(x, q, l, h, tab, prob, scheme))
    end, y)

    J = gy .+ gx * (-(Fx \ Fy))
    return J, qout, λout, nr, true
end

"defect in the symplecticity condition of one step"
function defect(prob, qₙ, λₙ, h, tab, scheme)
    J, qout, λout, nr, ok = stepmap(prob, qₙ, λₙ, h, tab, scheme)
    ok || return NaN, nr
    if extended(scheme)
        D = transpose(J) * Omega_lambda(prob, qout, λout, h) * J .- Omega_lambda(prob, qₙ, λₙ, h)
        return norm(D) / max(norm(Omega_lambda(prob, qₙ, λₙ, h)), eps()), nr
    else
        D = transpose(J) * Omega_bar(prob, qout) * J .- Omega_bar(prob, qₙ)
        return norm(D) / max(norm(Omega_bar(prob, qₙ)), eps()), nr
    end
end


# ---------------------------------------------------------------------------
# validation of the harness: the unprojected VPRK method is canonically symplectic
# ---------------------------------------------------------------------------

"""
One step of the unprojected VPRK method as a map on T*M, `(q,p) ↦ (q',p')`, for
initial data that need not satisfy the constraint.  Used only to validate the
differentiation machinery: the result must be canonically symplectic.
"""
function vprk_canonical(prob, z, h, tab)
    d = length(z) ÷ 2
    s = length(tab.b)
    res = function (V, zz)
        q, p = zz[1:d], zz[d+1:2d]
        Vs = [V[d*(i-1)+1:d*i] for i in 1:s]
        Q = [q .+ h .* sum(tab.a[i, j] .* Vs[j] for j in 1:s) for i in 1:s]
        F = [fϑ(prob, Q[i], Vs[i]) .- ∇H(prob, Q[i]) for i in 1:s]
        r = vcat([prob.ϑ(Q[i]) .- p .- h .* sum(tab.ā[i, j] .* F[j] for j in 1:s) for i in 1:s]...)
        znew = vcat(q .+ h .* sum(tab.b[i] .* Vs[i] for i in 1:s),
            p .+ h .* sum(tab.b[i] .* F[i] for i in 1:s))
        r, znew
    end
    V = repeat(qdot(prob, z[1:d]), s)
    nr = norm(res(V, z)[1])
    for _ in 1:100
        nr < 1e-14 && break
        Δ = ForwardDiff.jacobian(y -> res(y, z)[1], V) \ res(V, z)[1]
        all(isfinite, Δ) || break
        α = 1.0
        for _ in 1:30
            Vt = V .- α .* Δ
            nt = norm(res(Vt, z)[1])
            if isfinite(nt) && nt < nr
                V, nr = Vt, nt
                break
            end
            α /= 2
        end
    end
    isfinite(nr) && nr < 1e-10 || return NaN, nr
    Fx = ForwardDiff.jacobian(y -> res(y, z)[1], V)
    gx = ForwardDiff.jacobian(y -> res(y, z)[2], V)
    Fz = ForwardDiff.jacobian(y -> res(V, y)[1], z)
    gz = ForwardDiff.jacobian(y -> res(V, y)[2], z)
    J = gz .+ gx * (-(Fx \ Fz))
    Ω = [zeros(d, d) I(d); -I(d) zeros(d, d)]
    norm(transpose(J) * Ω * J .- Ω), norm(res(V, z)[1])
end


# ---------------------------------------------------------------------------
# the obstruction: a midpoint stage forces R(∞) = -1
# ---------------------------------------------------------------------------

"""
Verify the lemma behind the failure of the symplecticity proof of the midpoint
projection:  if some stage satisfies  Z_m = ½(z̄_n + z̄_{n+1}), i.e. row m of the
tableau equals bᵀ/2, then

    e_mᵀ A = ½ bᵀ  ⟹  bᵀA⁻¹ = 2 e_mᵀ  ⟹  bᵀA⁻¹e = 2  ⟹  R(∞) = 1 - bᵀA⁻¹e = -1 ,

so the midpoint-stage property and the R(∞) = +1 required by the proof are
mutually exclusive.  The same holds for the partner tableau Ā.
"""
function check_midpoint_stage_lemma()
    println("Midpoint-stage property vs. R(∞)")
    println("-"^78)
    @printf("%-8s %5s %10s %10s %14s %14s\n", "tableau", "s", "R(∞)[A]", "R(∞)[Ā]", "row m of A=b/2", "row m of Ā=b/2")
    for tab in tableaus()
        s = length(tab.b)
        R̄∞ = 1 - dot(tab.b, tab.ā \ ones(s))
        m = tab.midpointstage
        rowa = m > 0 ? "yes (m=$m)" : "no"
        rowā = m > 0 && isapprox(tab.ā[m, :], tab.b ./ 2; atol=1e-12) ? "yes (m=$m)" : "no"
        @printf("%-8s %5d %10.4f %10.4f %14s %14s\n", tab.name, s, tab.R∞, R̄∞, rowa, rowā)
        # symplecticity conditions of the tableau
        @assert all(isapprox(tab.b[i] * tab.ā[i, j] + tab.b[j] * tab.a[j, i], tab.b[i] * tab.b[j]; atol=1e-12)
                    for i in 1:s, j in 1:s)
        m > 0 && @assert isapprox(tab.R∞, -1.0; atol=1e-12) "midpoint stage must force R(∞) = -1"
    end
    println()
end


# ---------------------------------------------------------------------------
# reports
# ---------------------------------------------------------------------------

const STEPSIZES = [0.2, 0.1, 0.05, 0.025, 0.0125]

function orders(defects)
    o = fill(NaN, length(defects))
    for i in 2:length(defects)
        if isfinite(defects[i]) && isfinite(defects[i-1]) && defects[i] > 0
            o[i] = log2(defects[i-1] / defects[i])
        end
    end
    o
end

"""
Reference point for the Jacobian test: one step of the corresponding method from
the initial condition of the problem, so that q and λ lie on an actual
trajectory (λ = O(h) rather than an arbitrary value).
"""
function refpoint(prob, h, tab, scheme)
    qₙ, λₙ = copy(prob.q₀), zeros(length(prob.q₀))
    _, q, λ, _, ok = stepmap(prob, qₙ, λₙ, h, tab, scheme)
    ok ? (q, λ) : (qₙ, λₙ)
end

"""
    check_generalised_form(prob, tab, h)

The symplecticity analysis of the symmetric and symplectic projections rests on
the exact identity (`symmetric_projection_symplecticity_condition`)

    dp_n ∧ dq_n - dA_n ∧ dB_n  =  dp_{n+1} ∧ dq_{n+1} - dA_{n+1} ∧ dB_{n+1} ,
    A^i = h λ^k φ^k_{q^i}(z) = -h λ^k ϑ_{k,i}(q) ,   B^i = h λ^i ,

which for the symplectic projection (λ_n in the perturbation, λ_{n+1} in the
projection) is claimed to be the conservation of the modified form ω_λ on
M × R^d.  This routine evaluates *both* sides of that identity directly as
2-forms on the space of initial data (q_n, λ_n), so it distinguishes

  * the identity holds but ω_λ as printed is not the corresponding matrix
    (an algebra slip in the paper), from
  * the identity itself failing (the method is not symplectic in that sense).

Returns `(defect of the identity, defect of ω_λ preservation)`, both relative.
"""
function check_generalised_form(prob, tab, h)
    d = length(prob.q₀)
    qₙ, λₙ = refpoint(prob, h, tab, :symplectic)
    J, qₙ₊₁, λₙ₊₁, _, ok = stepmap(prob, qₙ, λₙ, h, tab, :symplectic)
    ok || return NaN, NaN

    # gradients (as row vectors on the 4-dimensional space of (q_n, λ_n))
    E = Matrix{Float64}(I, 2d, 2d)
    gq(i, which) = which === :n ? E[i, :] : J[i, :]
    gλ(k, which) = which === :n ? E[d+k, :] : J[d+k, :]

    function forms(q, λ, which)
        Jϑ = ForwardDiff.jacobian(prob.ϑ, q)                     # Jϑ[k,i] = ϑ_{k,i}
        Hϑ = [ForwardDiff.hessian(y -> prob.ϑ(y)[k], q) for k in 1:d]
        wedge(a, b) = a * transpose(b) .- b * transpose(a)
        dpdq = zeros(2d, 2d)
        for i in 1:d
            gp = sum(Jϑ[i, j] .* gq(j, which) for j in 1:d)       # dp^i = ϑ_{i,j} dq^j
            dpdq .+= wedge(gp, gq(i, which))
        end
        dAdB = zeros(2d, 2d)
        for i in 1:d
            gA = -h .* (sum(λ[k] .* Hϑ[k][i, j] .* gq(j, which) for k in 1:d, j in 1:d) .+
                        sum(Jϑ[k, i] .* gλ(k, which) for k in 1:d))
            gB = h .* gλ(i, which)
            dAdB .+= wedge(gA, gB)
        end
        dpdq .- dAdB
    end

    Sₙ = forms(qₙ, λₙ, :n)
    Sₙ₊₁ = forms(qₙ₊₁, λₙ₊₁, :np1)
    iddefect = norm(Sₙ₊₁ .- Sₙ) / max(norm(Sₙ), eps())

    # ω_λ with the corrected and with the printed sign of the mixed term
    defects = map((+1, -1)) do s
        Ωλ = Omega_lambda(prob, qₙ, λₙ, h; mixedsign=s)
        norm(transpose(J) * Omega_lambda(prob, qₙ₊₁, λₙ₊₁, h; mixedsign=s) * J .- Ωλ) /
        max(norm(Ωλ), eps())
    end

    iddefect, defects[1], defects[2]
end

function report_generalised_form(prob, tab)
    println("Symplectic projection: exact identity and the two signs of ω_λ (", tab.name, ")")
    println("-"^78)
    @printf("%-8s %20s %20s %20s\n", "h", "identity dp∧dq-dA∧dB", "ω_λ corrected sign", "ω_λ printed sign")
    for h in [0.2, 0.1, 0.05, 0.025]
        a, b, c = check_generalised_form(prob, tab, h)
        @printf("%-8s %20.3e %20.3e %20.3e\n", string(h), a, b, c)
    end
    println()
end

"""
Structural diagnostic explaining *why* a problem may fail to discriminate
between the projection methods.  For the constraint φ = p - ϑ(q) the defect
terms of all endpoint projections are (cf. the standard-projection section)

    (h²/2) Ω̄_ij dλ^i∧dλ^j   and   h² λ^k ϑ_{k,ij} dq^i∧dλ^j ,

so they vanish identically whenever

  * ϑ_k is affine for every component k with λ^k ≢ 0   (kills the second term), and
  * all but one component of λ vanish identically       (kills dλ^i∧dλ^j for d = 2).

Both Lotka-Volterra models satisfy this: their drift is confined to the first
constraint component, so λ = (0, λ²), while ϑ₂ is affine (ϑ₂ = q₁ resp. ϑ₂ ≡ 0).
Every projection method is then exactly symplectic on these problems -- including
the standard projection, which is provably not symplectic in general.  Such
problems therefore cannot be used to test symplecticity.
"""
function report_degeneracy(prob, tab; h=0.1)
    d = length(prob.q₀)
    x, nr, ok = newton(prob, prob.q₀, zeros(d), h, tab, :standard)
    _, _, λ = residual(x, prob.q₀, zeros(d), h, tab, prob, :standard)
    hess = [norm(ForwardDiff.hessian(y -> prob.ϑ(y)[k], prob.q₀)) for k in 1:d]
    active = [abs(λ[k]) > 1e-10 for k in 1:d]
    kills2 = all(!active[k] || hess[k] < 1e-12 for k in 1:d)
    kills1 = count(active) <= 1
    @printf("%-26s λ = %-26s ‖∂²ϑ_k‖ = %-22s %s\n", prob.name,
        string(round.(λ; sigdigits=3)), string(round.(hess; sigdigits=3)),
        (kills1 && kills2) ? "DEGENERATE (defect vanishes structurally)" : "discriminating")
    ok || @printf("%-26s (warning: step did not converge, ‖res‖ = %.2e)\n", "", nr)
end

function report_defects(prob, tab)
    @printf("%s / %s   (R(∞) = %+.0f%s)\n", prob.name, tab.name, tab.R∞,
        tab.midpointstage > 0 ? ", midpoint stage m=$(tab.midpointstage)" : "")
    println("-"^78)
    @printf("%-13s %-9s", "scheme", "form")
    for h in STEPSIZES
        @printf("%11s", "h=$h")
    end
    @printf("   %s\n", "order")
    for scheme in ALL_SCHEMES
        defects = Float64[]
        for h in STEPSIZES
            q, λ = refpoint(prob, h, tab, scheme)
            d, _ = defect(prob, q, λ, h, tab, scheme)
            push!(defects, d)
        end
        o = orders(defects)
        @printf("%-13s %-9s", string(scheme), extended(scheme) ? "ω_λ" : "ω̄")
        for d in defects
            isfinite(d) ? @printf("%11.2e", d) : @printf("%11s", "---")
        end
        oo = filter(isfinite, o)
        @printf("   %s\n", isempty(oo) ? "unsolvable" : @sprintf("%.2f", last(oo)))
    end
    println()
end

"""
Reference solution of the reduced Euler-Lagrange equations Ω̄(q) q̇ = ∇H(q) by
classical Runge-Kutta.  Also validates that the projected schemes converge to
the *correct* dynamics, not merely to something.
"""
function reference(prob, T; n=200000)
    q = copy(prob.q₀)
    dt = T / n
    for _ in 1:n
        k1 = qdot(prob, q)
        k2 = qdot(prob, q .+ (dt / 2) .* k1)
        k3 = qdot(prob, q .+ (dt / 2) .* k2)
        k4 = qdot(prob, q .+ dt .* k3)
        q = q .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    q
end

function report_convergence(prob, tab)
    println("Convergence to the exact flow (global error at T = 1, ", tab.name, ")")
    println("-"^78)
    qref = reference(prob, 1.0)
    @printf("%-13s", "scheme")
    for h in [0.1, 0.05, 0.025]
        @printf("%13s", "h=$h")
    end
    @printf("   %s\n", "order")
    for scheme in ALL_SCHEMES
        errs = Float64[]
        for h in [0.1, 0.05, 0.025]
            q, λ = copy(prob.q₀), zeros(length(prob.q₀))
            good = true
            for _ in 1:round(Int, 1 / h)
                _, q, λ, _, ok = stepmap(prob, q, λ, h, tab, scheme)
                if !ok
                    good = false
                    break
                end
            end
            push!(errs, good ? norm(q .- qref) : NaN)
        end
        o = orders(errs)
        @printf("%-13s", string(scheme))
        for e in errs
            isfinite(e) ? @printf("%13.3e", e) : @printf("%13s", "---")
        end
        oo = filter(isfinite, o)
        @printf("   %s\n", isempty(oo) ? "unsolvable" : @sprintf("%.2f", last(oo)))
    end
    println()
end

"""
Demonstrate that the R(∞) = +1 midpoint projection (the sign required by the
symplecticity proof) has no solution: with that sign the perturbation cancels
out of the midpoint,  z̄_{n+1/2} = ½(z_n + z_{n+1}),  so the constraint that
holds automatically at the internal stage degenerates into

    ½ (ϑ(q_n) + ϑ(q_{n+1})) = ϑ(½ (q_n + q_{n+1})) ,

a condition that involves neither λ nor H and is generically violated.  The
irreducible residual of the step equations equals 2κ with
κ = ϑ(q_mid) - ½(ϑ(q_n) + ϑ(q_{n+1})).
"""
function report_midpoint_R1_obstruction(prob, tab; ntrials=60)
    println("Solvability of the midpoint projection with the R(∞) = +1 sign")
    println("-"^78)
    @printf("%-8s %14s %14s %14s\n", "h", "best ‖res‖", "|2κ| there", "‖res‖ (R(∞)=-1)")
    d = length(prob.q₀)
    for h in [0.2, 0.1, 0.05]
        best, bq = Inf, nothing
        for t in 1:ntrials
            x = initialguess(prob, prob.q₀, h, tab, :midpoint_R1) .+ (t == 1 ? 0.0 : 0.3) .* randn(nunknowns(:midpoint_R1, length(tab.b), d))
            f = y -> residual(y, prob.q₀, zeros(d), h, tab, prob, :midpoint_R1)[1]
            nr = norm(f(x))
            for _ in 1:80
                J = ForwardDiff.jacobian(f, x)
                Δ = try
                    J \ f(x)
                catch
                    break
                end
                all(isfinite, Δ) || break
                α = 1.0
                improved = false
                for _ in 1:25
                    xt = x .- α .* Δ
                    nt = norm(f(xt))
                    if isfinite(nt) && nt < nr
                        x, nr, improved = xt, nt, true
                        break
                    end
                    α /= 2
                end
                improved || break
            end
            if nr < best
                best = nr
                bq = residual(x, prob.q₀, zeros(d), h, tab, prob, :midpoint_R1)[2]
            end
        end
        κ = prob.ϑ((prob.q₀ .+ bq) ./ 2) .- (prob.ϑ(prob.q₀) .+ prob.ϑ(bq)) ./ 2
        _, nrm, _ = newton(prob, prob.q₀, zeros(d), h, tab, :midpoint)
        @printf("%-8s %14.3e %14.3e %14.3e\n", string(h), best, norm(2 .* κ), nrm)
    end
    println()
end


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

function main()
    println("="^78)
    println("Symplecticity of projected variational integrators for degenerate")
    println("Lagrangian systems -- numerical verification")
    println("="^78)
    println()

    check_midpoint_stage_lemma()

    println("Validation of the harness: unprojected VPRK on T*M (off the constraint)")
    println("-"^78)
    @printf("%-26s %-8s %14s %14s\n", "problem", "tableau", "|JᵀΩJ-Ω|", "‖res‖")
    for prob in testproblems(), tab in tableaus()
        z = vcat(prob.q₀, prob.ϑ(prob.q₀) .+ [0.05, -0.03])   # deliberately off Δ
        dfc, nr = vprk_canonical(prob, z, 0.05, tab)
        @printf("%-26s %-8s %14.3e %14.3e\n", prob.name, tab.name, dfc, nr)
    end
    println()

    println("Is the problem able to discriminate between the projection methods?")
    println("-"^78)
    for prob in testproblems()
        report_degeneracy(prob, tableaus()[1])
    end
    println()

    for prob in testproblems()
        println("#"^78)
        println("# ", prob.name)
        println("#"^78)
        println()
        for tab in tableaus()
            report_defects(prob, tab)
        end
        report_convergence(prob, tableaus()[1])
        report_generalised_form(prob, tableaus()[1])
        report_midpoint_R1_obstruction(prob, tableaus()[1])
    end
end

end # module

isinteractive() || ProjectedVPRKSymplecticity.main()

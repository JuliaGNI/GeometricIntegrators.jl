@doc raw"""
Discontinuous Galerkin Variational Integrator with a projection at the **initial** time
of each interval. *EXPERIMENTAL.*

The interior contribution and the trapezoidal flux are those of [`DGVI`](@ref), but the
step is closed by a projection condition imposed at ``t_n`` rather than by the
definition of ``p_n``,
```math
\vartheta (q_n^+) - \vartheta (q_n^-) - \nabla \vartheta^{T} (q_n) \cdot (q_n^+ - q_n^-) = 0 ,
```
in which the one-sided one-forms are reconstructed from the *stage* values of the
one-form,
```math
\vartheta (q_n^+) = \sum \limits_{i} r^{+}_{i} P_{n,i} ,
\qquad
\vartheta (q_{n+1}^-) = \sum \limits_{i} r^{-}_{i} P_{n,i} .
```
Consequently ``\vartheta (q_n^-)`` and ``q_n^-`` are carried over from the previous step.

!!! warning
    That reconstruction applies the *position* basis' boundary coefficients to values
    sampled at the *quadrature* nodes, which is only meaningful when the two node sets
    coincide. The constructor therefore requires `nbasis(basis) == nnodes(quadrature)`,
    as holds for `Lagrange(QuadratureRules.nodes(quadrature))`. The corresponding legacy
    code looped `for i in 1:S` over an array of length `R` and would raise a
    `BoundsError` whenever `S > R`.

!!! note
    This variant has no written derivation in the manuscript, and its momentum output is
    diagnostic: the legacy implementation never set `p`, and here
    ``p_{n+1} = \vartheta (q_{n+1})`` is supplied so that the solution is complete.

### Constructor

```
DGVIP0(basis::Basis, quadrature::QuadratureRule)
```
"""
DGVIP0

GeometricBase.description(::DGVIP0) = "Discontinuous Galerkin Variational Integrator with Initial Projection"


function jump!(sol, params, int::GeometricIntegrator{<:DGVIP0}, ST)
    local C = cache(int, ST)
    local M = method(int)
    local st = cache(int).state
    local t₀ = sol.t - timestep(int)
    local t₁ = sol.t

    # jumps, using the value carried over from the previous step
    C.λ .= C.q⁺ .- st.q⁻
    C.λ⁺ .= C.q⁺ .- sol.q
    C.λ̄⁻ .= C.q̄ .- C.q̄⁻

    equations(int).ϑ(C.θ, t₀, sol.q, sol.q, params)
    equations(int).ϑ(C.Θ̅, t₁, C.q̄, C.q̄, params)

    # boundary one-forms reconstructed from the stage values
    for k in eachindex(C.θ⁺, C.Θ̅⁻)
        y⁺ = y⁻ = zero(ST)
        for i in eachindex(C.P)
            y⁺ += M.r⁺[i] * C.P[i][k]
            y⁻ += M.r⁻[i] * C.P[i][k]
        end
        C.θ⁺[k] = y⁺
        C.Θ̅⁻[k] = y⁻
    end

    equations(int).g(C.g, t₀, sol.q, C.z, C.λ, params)
    equations(int).g(C.g⁺, t₀, C.q⁺, C.z, C.λ⁺, params)
    equations(int).g(C.ḡ⁻, t₁, C.q̄⁻, C.z, C.λ̄⁻, params)

    C.p̄ .= C.Θ̅
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DGVIP0}) where {ST}
    local C = cache(int, ST)
    local M = method(int)
    local D = length(C.q̃)
    local S = nbasis(M)

    residual_interior!(b, sol, params, int)

    for i in eachindex(M.r⁻, M.r⁺)
        for k in 1:D
            b[D*(i-1)+k] += M.r⁺[i] * (C.θ[k] + C.θ⁺[k]) / 2
            b[D*(i-1)+k] -= M.r⁻[i] * (C.Θ̅[k] + C.Θ̅⁻[k]) / 2
            b[D*(i-1)+k] += M.r⁺[i] * C.g⁺[k] / 2
            b[D*(i-1)+k] += M.r⁻[i] * C.ḡ⁻[k] / 2
        end
    end

    # closure: the projection at the initial time
    for k in 1:D
        b[D*S+k] = C.θ⁺[k] - cache(int).state.θ⁻[k] - C.g[k]
    end
end


function update_state!(int::GeometricIntegrator{<:DGVIP0}, DT)
    local C = cache(int, DT)
    local st = cache(int).state
    st.q⁻ .= C.q̄⁻
    st.θ⁻ .= C.Θ̅⁻
    return
end

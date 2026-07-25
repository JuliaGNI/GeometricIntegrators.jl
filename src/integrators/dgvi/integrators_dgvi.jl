@doc raw"""
Discontinuous Galerkin Variational Integrator with a trapezoidal numerical flux.

The discrete action augments the interior quadrature with the flux evaluated as a
trapezoidal rule between the nodal value ``q_n`` and the two one-sided limits,
```math
\mathcal{A}_d = h \sum \limits_{n} \bigg[ \sum \limits_{i=1}^{R} b_i \big[ \vartheta (Q_{n,i}) \cdot V_{n,i} - H(Q_{n,i}) \big]
 + \frac{\vartheta (q_n) + \vartheta (q_n^+)}{2} \cdot (q_n^+ - q_n)
 + \frac{\vartheta (q_{n+1}^-) + \vartheta (q_{n+1})}{2} \cdot (q_{n+1} - q_{n+1}^-) \bigg] ,
```
whose stationarity gives, for all ``j``,
```math
0 = \sum \limits_{i=1}^{R} b_i \big[ h \, m_{ij} \, F_{n,i} + a_{ij} \, P_{n,i} \big]
 + r^{+}_{j} \, \frac{\vartheta (q_n) + \vartheta (q_n^+)}{2}
 - r^{-}_{j} \, \frac{\vartheta (q_{n+1}^-) + \vartheta (q_{n+1})}{2}
 + \frac{r^{+}_{j}}{2} \, \nabla \vartheta^{T} (q_n^+) \cdot (q_n^+ - q_n)
 + \frac{r^{-}_{j}}{2} \, \nabla \vartheta^{T} (q_{n+1}^-) \cdot (q_{n+1} - q_{n+1}^-) ,
```
together with the variation with respect to the nodal value,
```math
\vartheta (q_n^+) = \vartheta (q_n^-) + \nabla \vartheta^{T} (q_n) \cdot (q_n^+ - q_n^-) .
```
Introducing
```math
p_{n} = \vartheta (q_{n}^-) + \nabla \vartheta^{T} (q_{n}) \cdot (q_{n} - q_{n}^-)
```
turns the latter into
```math
\vartheta (q_n^+) = p_n + \nabla \vartheta^{T} (q_n) \cdot (q_n^+ - q_n) ,
```
so the scheme is a genuine map ``(q_n, p_n) \mapsto (q_{n+1}, p_{n+1})`` and needs no
state beyond the solution itself — unlike the other four variants. The step is closed by
```math
p_{n+1} = \vartheta (q_{n+1}^-) + \nabla \vartheta^{T} (q_{n+1}) \cdot (q_{n+1} - q_{n+1}^-) .
```
Initial conditions ``q_0`` and ``p_0 = \vartheta(q_0)`` must be prescribed.

The unknowns are the ``S`` basis coefficients plus ``q_{n+1}``.

### Constructor

```
DGVI(basis::Basis, quadrature::QuadratureRule)
```
"""
DGVI

GeometricBase.description(::DGVI) = "Discontinuous Galerkin Variational Integrator"

# `DGVI` is the one variant that is a genuine (q,p) map: `p` carries the jump
# information, so nothing has to be propagated separately.
initialise_state!(sol, params, ::GeometricIntegrator{<:DGVI}) = nothing


function jump!(sol, params, int::GeometricIntegrator{<:DGVI}, ST)
    local C = cache(int, ST)
    local t₀ = sol.t - timestep(int)
    local t₁ = sol.t

    # jumps at the two ends of the interval
    C.λ⁺ .= C.q⁺ .- sol.q
    C.λ̄⁻ .= C.q̄ .- C.q̄⁻

    # one-form at the nodal values and the one-sided limits
    equations(int).ϑ(C.θ, t₀, sol.q, sol.q, params)
    equations(int).ϑ(C.θ⁺, t₀, C.q⁺, C.q⁺, params)
    equations(int).ϑ(C.Θ̅, t₁, C.q̄, C.q̄, params)
    equations(int).ϑ(C.Θ̅⁻, t₁, C.q̄⁻, C.q̄⁻, params)

    # ∇ϑᵀ contracted with the jumps
    equations(int).g(C.g⁺, t₀, C.q⁺, C.z, C.λ⁺, params)
    equations(int).g(C.ḡ⁻, t₁, C.q̄⁻, C.z, C.λ̄⁻, params)
    equations(int).g(C.h⁺, t₀, sol.q, C.z, C.λ⁺, params)
    equations(int).g(C.h̅⁻, t₁, C.q̄, C.z, C.λ̄⁻, params)

    # pₙ₊₁ = ϑ(qₙ₊₁⁻) + ∇ϑᵀ(qₙ₊₁)⋅(qₙ₊₁ - qₙ₊₁⁻)
    C.p̄ .= C.Θ̅⁻ .+ C.h̅⁻
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DGVI}) where {ST}
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

    # closure: ϑ(qₙ⁺) - pₙ - ∇ϑᵀ(qₙ)⋅(qₙ⁺ - qₙ) = 0
    for k in 1:D
        b[D*S+k] = C.θ⁺[k] - sol.p[k] - C.h⁺[k]
    end
end

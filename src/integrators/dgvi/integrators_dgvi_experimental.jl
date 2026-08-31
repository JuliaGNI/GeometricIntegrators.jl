@doc raw"""
Discontinuous Galerkin Variational Integrator with a **midpoint-average** numerical
flux. *EXPERIMENTAL.*

Instead of the trapezoidal flux of [`DGVI`](@ref), the one-form is evaluated once, at
the average of the two one-sided limits,
```math
\langle q \rangle_{n} = \tfrac{1}{2} \big( q_{n}^{-} + q_{n}^{+} \big) ,
\qquad
[\![ q ]\!]_{n} = q_{n}^{+} - q_{n}^{-} ,
```
giving the flux ``\vartheta ( \langle q \rangle_{n} ) \cdot [\![ q ]\!]_{n}`` — the
jump/average form of the manuscript with averaging parameter ``\alpha = 1/2``. The
discrete equations become
```math
0 = \sum \limits_{i=1}^{R} b_i \big[ h \, m_{ij} \, F_{n,i} + a_{ij} \, P_{n,i} \big]
 + r^{+}_{j} \, \vartheta ( \langle q \rangle_{n} )
 - r^{-}_{j} \, \vartheta ( \langle q \rangle_{n+1} )
 + \tfrac{1}{2} r^{+}_{j} \, \nabla \vartheta^{T} ( \langle q \rangle_{n} ) \cdot [\![ q ]\!]_{n}
 + \tfrac{1}{2} r^{-}_{j} \, \nabla \vartheta^{T} ( \langle q \rangle_{n+1} ) \cdot [\![ q ]\!]_{n+1} ,
```
closed by identifying the nodal value with the average,
```math
q_{n+1} = \langle q \rangle_{n+1} ,
```
and by continuity with the carried-over right limit, ``q_{n}^{+} - r^{+} \cdot x_{n} = 0``.
Both ``q_{n}^{-}`` and ``q_{n}^{+}`` are propagated from step to step.

!!! note
    The most heuristic of the five variants: no derivation exists for it, and its
    momentum output is diagnostic — ``p_{n+1} = \vartheta ( \langle q \rangle_{n+1} )``.

### Constructor

```
DGVIEXP(basis::Basis, quadrature::QuadratureRule)
```
"""
DGVIEXP

function GeometricBase.description(::DGVIEXP)
    "Discontinuous Galerkin Variational Integrator (midpoint-average flux)"
end

# two trailing blocks: q_{n+1} and q_{n+1}⁺
nclosure(::DGVIEXP) = 2

function jump!(sol, params, int::GeometricIntegrator{<:DGVIEXP}, ST)
    local C = cache(int, ST)
    local st = dgvi_state(int)
    local t₀ = sol.t - timestep(int)
    local t₁ = sol.t

    # averages and jumps across the two interval boundaries
    C.ϕ .= (st.q⁻ .+ C.q⁺) ./ 2
    C.ϕ̅ .= (C.q̄⁻ .+ C.q̄⁺) ./ 2
    C.λ .= C.q⁺ .- st.q⁻
    C.λ̄ .= C.q̄⁺ .- C.q̄⁻

    equations(int).ϑ(C.θ, t₀, C.ϕ, C.ϕ, params)
    equations(int).ϑ(C.Θ̅, t₁, C.ϕ̅, C.ϕ̅, params)

    equations(int).g(C.g, t₀, C.ϕ, C.z, C.λ, params)
    equations(int).g(C.ḡ, t₁, C.ϕ̅, C.z, C.λ̄, params)

    C.p̄ .= C.Θ̅
end

function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DGVIEXP}) where {ST}
    local C = cache(int, ST)
    local M = method(int)
    local D = length(C.q̃)
    local S = nbasis(M)
    local st = dgvi_state(int)

    residual_interior!(b, sol, params, int)

    for i in eachindex(M.r⁻, M.r⁺)
        for k in 1:D
            b[D * (i - 1) + k] += M.r⁺[i] * C.θ[k]
            b[D * (i - 1) + k] -= M.r⁻[i] * C.Θ̅[k]
            b[D * (i - 1) + k] += M.r⁺[i] * C.g[k] / 2
            b[D * (i - 1) + k] += M.r⁻[i] * C.ḡ[k] / 2
        end
    end

    # the nodal value is the average of the one-sided limits
    for k in 1:D
        b[D * S + k] = C.q̄[k] - C.ϕ̅[k]
    end

    # continuity with the carried-over right limit
    for k in 1:D
        b[D * (S + 1) + k] = st.q⁺[k] - C.q⁺[k]
    end
end

function update_state!(int::GeometricIntegrator{<:DGVIEXP}, DT)
    local C = cache(int, DT)
    local st = dgvi_state(int)
    st.q⁻ .= C.q̄⁻
    st.q⁺ .= C.q̄⁺
    return
end

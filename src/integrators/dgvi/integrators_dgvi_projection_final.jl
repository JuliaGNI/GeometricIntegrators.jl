@doc raw"""
Discontinuous Galerkin Variational Integrator with a projection at the **final** time of
each interval. *EXPERIMENTAL.*

The interior contribution and the trapezoidal flux are those of [`DGVI`](@ref). The step
is closed by two conditions: continuity of the carried-over right limit,
```math
q_{n}^{+} - r^{+} \cdot x_{n} = 0 ,
```
and a projection imposed at ``t_{n+1}``,
```math
\vartheta (q_{n+1}^{+}) - \vartheta (q_{n+1}^{-}) - \nabla \vartheta^{T} (q_{n+1}) \cdot (q_{n+1}^{+} - q_{n+1}^{-}) = 0 .
```
The unknowns are therefore the ``S`` basis coefficients, ``q_{n+1}`` and
``q_{n+1}^{+}``, and ``q_{n}^{+}`` is carried over from the previous step.

!!! note
    No written derivation exists for this variant. The manuscript's Appendix A flags
    precisely this family as *"under-determined by 2d equations"*, closed only by
    postulating a relation between the nodal value and the one-sided limits, with the
    author's own footnote *"This story is not really clean"*. Its momentum output is
    diagnostic: ``p_{n+1} = \vartheta (q_{n+1}^{+})``.

### Constructor

```
DGVIP1(basis::Basis, quadrature::QuadratureRule)
```
"""
DGVIP1

GeometricBase.description(::DGVIP1) = "Discontinuous Galerkin Variational Integrator with Final Projection"

# two trailing blocks: q_{n+1} and q_{n+1}⁺
nclosure(::DGVIP1) = 2


function jump!(sol, params, int::GeometricIntegrator{<:DGVIP1}, ST)
    local C = cache(int, ST)
    local t₀ = sol.t - timestep(int)
    local t₁ = sol.t

    C.λ⁺ .= C.q⁺ .- sol.q
    C.λ̄⁻ .= C.q̄ .- C.q̄⁻
    C.λ̄ .= C.q̄⁺ .- C.q̄⁻

    equations(int).ϑ(C.θ, t₀, sol.q, sol.q, params)
    equations(int).ϑ(C.θ⁺, t₀, C.q⁺, C.q⁺, params)
    equations(int).ϑ(C.Θ̅, t₁, C.q̄, C.q̄, params)
    equations(int).ϑ(C.Θ̅⁻, t₁, C.q̄⁻, C.q̄⁻, params)
    equations(int).ϑ(C.Θ̅⁺, t₁, C.q̄⁺, C.q̄⁺, params)

    equations(int).g(C.g⁺, t₀, C.q⁺, C.z, C.λ⁺, params)
    equations(int).g(C.ḡ⁻, t₁, C.q̄⁻, C.z, C.λ̄⁻, params)
    equations(int).g(C.ḡ, t₁, C.q̄, C.z, C.λ̄, params)

    C.p̄ .= C.Θ̅⁺
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DGVIP1}) where {ST}
    local C = cache(int, ST)
    local M = method(int)
    local D = length(C.q̃)
    local S = nbasis(M)
    local st = dgvi_state(int)

    residual_interior!(b, sol, params, int)

    for i in eachindex(M.r⁻, M.r⁺)
        for k in 1:D
            b[D*(i-1)+k] += M.r⁺[i] * (C.θ[k] + C.θ⁺[k]) / 2
            b[D*(i-1)+k] -= M.r⁻[i] * (C.Θ̅[k] + C.Θ̅⁻[k]) / 2
            b[D*(i-1)+k] += M.r⁺[i] * C.g⁺[k] / 2
            b[D*(i-1)+k] += M.r⁻[i] * C.ḡ⁻[k] / 2
        end
    end

    # continuity with the carried-over right limit
    for k in 1:D
        b[D*S+k] = st.q⁺[k] - C.q⁺[k]
    end

    # projection at the final time
    for k in 1:D
        b[D*(S+1)+k] = C.Θ̅⁺[k] - C.Θ̅⁻[k] - C.ḡ[k]
    end
end


function update_state!(int::GeometricIntegrator{<:DGVIP1}, DT)
    dgvi_state(int).q⁺ .= cache(int, DT).q̄⁺
    return
end

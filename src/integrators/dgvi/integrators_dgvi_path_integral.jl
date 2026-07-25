@doc raw"""
Discontinuous Galerkin Variational Integrator with a **path-integral** numerical flux.

This is the general variant, and the one derived in `docs/src/integrators/dgvi.md`. The
non-conservative product of the one-form and the jump is discretised by integrating
along a path ``\Phi(\tau; q^-, q^+)`` connecting the two one-sided limits,
```math
\int \limits_0^1 \vartheta \big( \Phi(\tau; q^-, q^+) \big) \cdot \frac{d \Phi}{d\tau} \, d\tau
\approx \sum \limits_{i=1}^{F} \beta_i \, \vartheta \big( \Phi (\gamma_i) \big) \cdot \Phi' (\gamma_i) ,
```
with a quadrature rule of ``F`` nodes ``\gamma_i`` and weights ``\beta_i``. Writing the
path as ``\Phi(\tau) = q^- \varphi^-(\tau) + q^+ \varphi^+(\tau)`` and setting
```math
\mu^{\mp}_i = \varphi^{\mp} (\gamma_i) ,
\qquad
\alpha^{\mp}_i = \frac{d\varphi^{\mp}}{d\tau} (\gamma_i) ,
```
the path values and their derivatives are
```math
\Phi_i = \mu^-_i q^- + \mu^+_i q^+ ,
\qquad
\Lambda_i = \alpha^-_i q^- + \alpha^+_i q^+ ,
```
and the discrete Euler-Lagrange equations read
```math
0 = \sum \limits_{i=1}^{R} b_i \big[ h \, m_{ij} \, F_{n,i} + a_{ij} \, P_{n,i} \big]
 + \sum \limits_{i=1}^{F} \beta_i \Big[
     r^{+}_{j} \alpha^{+}_i \, \vartheta (\Phi_i)
   + r^{-}_{j} \alpha^{-}_i \, \vartheta (\bar{\Phi}_i)
   + r^{+}_{j} \mu^{+}_i \, \nabla \vartheta^{T} (\Phi_i) \cdot \Lambda_i
   + r^{-}_{j} \mu^{-}_i \, \nabla \vartheta^{T} (\bar{\Phi}_i) \cdot \bar{\Lambda}_i
 \Big] ,
```
where barred quantities are those of the jump at ``t_{n+1}``. Note that the full path is
integrated at *both* boundaries, but only the ``\partial / \partial q_n^+`` terms are
taken at ``t_n`` and only the ``\partial / \partial q_{n+1}^-`` terms at ``t_{n+1}``, so
there is no double counting.

The step is closed by reconstructing the nodal value at the midpoint of the path,
```math
q_{n} = \Phi (1/2; q_n^-, q_n^+) = \rho^- q_n^- + \rho^+ q_n^+ ,
```
and ``q_n^-`` is carried over from the previous step. With
`Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))` this reduces to the
trapezoidal flux on the full jump with ``\rho^- = \rho^+ = 1/2``.

!!! note
    The gauge term ``\nu`` of `docs/src/integrators/dgvi.md` is *not* implemented: the
    ``\rho^{\mp}`` above reconstruct the central value and are unrelated to it. Adding it
    would introduce two further flux terms; `GeometricProblems.LotkaVolterra2dGauge`
    exists to exercise a gauge-transformed problem in the meantime.

    Its momentum output is diagnostic: ``p_{n+1} = \vartheta (q_{n+1})``.

### Constructor

```
DGVIPI(basis::Basis, quadrature::QuadratureRule, jump::Discontinuity)
```
"""
DGVIPI

GeometricBase.description(::DGVIPI) = "Discontinuous Galerkin Variational Integrator (path-integral flux)"


"""
Seed the carried-over jump values. `q₀⁻ = q₀` and `q₀⁺` follows from the midpoint
reconstruction `q₀ = ρ⁻ q₀⁻ + ρ⁺ q₀⁺`.
"""
function initialise_state!(sol, params, int::GeometricIntegrator{<:DGVIPI})
    local st = cache(int).state
    local M = method(int)
    st.initialised && return
    st.q⁻ .= sol.q
    st.q⁺ .= (sol.q .- M.ρ⁻ .* st.q⁻) ./ M.ρ⁺
    st.θ⁻ .= sol.p
    st.initialised = true
    return
end


function jump!(sol, params, int::GeometricIntegrator{<:DGVIPI}, ST)
    local C = cache(int, ST)
    local M = method(int)
    local st = cache(int).state
    local t₀ = sol.t - timestep(int)
    local t₁ = sol.t

    # the right limit at t_{n+1} follows from the midpoint reconstruction
    C.q̄⁺ .= (C.q̄ .- M.ρ⁻ .* C.q̄⁻) ./ M.ρ⁺

    # path values and derivatives at the two boundaries
    for i in eachindex(C.Φ, C.Λ, C.Φ̄, C.Λ̄)
        C.Φ[i] .= M.μ⁻[i] .* st.q⁻ .+ M.μ⁺[i] .* C.q⁺
        C.Λ[i] .= M.α⁻[i] .* st.q⁻ .+ M.α⁺[i] .* C.q⁺
        C.Φ̄[i] .= M.μ⁻[i] .* C.q̄⁻ .+ M.μ⁺[i] .* C.q̄⁺
        C.Λ̄[i] .= M.α⁻[i] .* C.q̄⁻ .+ M.α⁺[i] .* C.q̄⁺
    end

    for i in eachindex(C.Θ, C.Θ̄, C.G, C.Ḡ)
        equations(int).ϑ(C.Θ[i], t₀, C.Φ[i], C.Λ[i], params)
        equations(int).ϑ(C.Θ̄[i], t₁, C.Φ̄[i], C.Λ̄[i], params)
        equations(int).g(C.G[i], t₀, C.Φ[i], C.z, C.Λ[i], params)
        equations(int).g(C.Ḡ[i], t₁, C.Φ̄[i], C.z, C.Λ̄[i], params)
    end

    equations(int).ϑ(C.p̄, t₁, C.q̄, C.q̄, params)
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DGVIPI}) where {ST}
    local C = cache(int, ST)
    local M = method(int)
    local D = length(C.q̃)
    local S = nbasis(M)
    local st = cache(int).state

    residual_interior!(b, sol, params, int)

    for i in eachindex(M.r⁻, M.r⁺)
        for k in 1:D
            z = zero(ST)
            for j in eachindex(C.Θ, C.Θ̄, C.G, C.Ḡ)
                z += M.β[j] * M.r⁺[i] * M.α⁺[j] * C.Θ[j][k]
                z += M.β[j] * M.r⁻[i] * M.α⁻[j] * C.Θ̄[j][k]
                z += M.β[j] * M.r⁺[i] * M.μ⁺[j] * C.G[j][k]
                z += M.β[j] * M.r⁻[i] * M.μ⁻[j] * C.Ḡ[j][k]
            end
            b[D*(i-1)+k] += z
        end
    end

    # closure: qₙ = Φ(1/2; qₙ⁻, qₙ⁺)
    for k in 1:D
        b[D*S+k] = sol.q[k] - M.ρ⁻ * st.q⁻[k] - M.ρ⁺ * C.q⁺[k]
    end
end


function update_state!(int::GeometricIntegrator{<:DGVIPI}, DT)
    cache(int).state.q⁻ .= cache(int, DT).q̄⁻
    return
end

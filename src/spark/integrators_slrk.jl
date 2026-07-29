"Holds all parameters of an Specialised Partitioned Additive Runge-Kutta method for variational systems subject to constraints."
struct SLRK{DT<:Number,DVT} <: LSPARKMethod
    name::Symbol
    o::Int
    s::Int
    r::Int

    q::Tableau{DT}
    p::Tableau{DT}

    q̃::Tableau{DT}
    p̃::Tableau{DT}

    ω::Matrix{DT}
    d::DVT

    function SLRK(name::Symbol, o::Int, s::Int,
        q::Tableau{DT}, p::Tableau{DT},
        q̃::Tableau{DT}, p̃::Tableau{DT},
        ω::Matrix{DT}, d::DVT=nothing) where {DT,DVT<:Union{AbstractVector,Nothing}}

        @assert s > 0 "Number of stages s must be > 0"

        @assert s == q.s == p.s == q̃.s == p̃.s
        @assert size(ω, 1) == s
        @assert size(ω, 2) == s + 1

        @assert d === nothing || length(d) == s

        new{DT,DVT}(name, o, s, s, q, p, q̃, p̃, ω, d)
    end
end

tableau(method::SLRK) = method

nstages(method::SLRK) = method.s
pstages(method::SLRK) = method.r

hasnullvector(method::SLRK{DT,Nothing}) where {DT} = false
hasnullvector(method::SLRK{DT,<:AbstractVector}) where {DT} = true

solversize(problem::LDAEProblem, method::SLRK) =
    4 * length(vec(initial_conditions(problem).q)) * nstages(method)


@doc raw"""
Specialised Lobatto Runge-Kutta integrator for degenerate variational systems
with projection on the secondary constraint.

`SLRK` targets degenerate Lagrangians of the form
```math
L (q, v) = \vartheta (q) \cdot v - H (q) ,
```
whose fibre derivative ``p = \vartheta(q)`` is not invertible for ``v``, so that the
equations of motion take the form of an index-two DAE with the *primary* constraint
``\phi (q,p) = p - \vartheta (q) = 0`` and the *secondary* constraint
``\psi = \dot{p} - \dot{q} \cdot \nabla \vartheta (q) = \tfrac{d}{dt} \phi = 0``.

Unlike the SPARK methods of Jay, which average the primary constraints over the
internal stages, `SLRK` imposes the primary constraint at *every* stage and averages
the *secondary* constraints instead. This is what makes the primary constraint hold
to round-off along the whole trajectory.

In contrast to [`VSPARKsecondary`](@ref) there is only a **single set of ``s`` stages**:
the projective stages coincide with the internal stages, so ``\sigma = s`` and the
``\Lambda`` multipliers live on the same nodes as the velocities ``V``.

One step ``(q_n, p_n) \mapsto (q_{n+1}, p_{n+1})`` solves
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a^1_{ij} V_{n,j}
                 + h \sum \limits_{j=1}^{s} a^2_{ij} \Lambda_{n,j} , & i &= 1, ..., s , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \bar{a}^1_{ij} F_{n,j}
                 + h \sum \limits_{j=1}^{s} \bar{a}^2_{ij} G_{n,j} , & i &= 1, ..., s , \\
0 &= \Phi_{n,i} - \frac{d_i}{\bar{b}^1_i} \, \mu , & i &= 1, ..., s , \\
0 &= \sum \limits_{j=1}^{s} \omega_{ij} \, \Psi_{n,j} + \omega_{i,s+1} \, \phi (q_{n+1}, p_{n+1}) ,
  & i &= 1, ..., s , \\
0 &= \sum \limits_{i=1}^{s} d_{i} \, V_{n,i} ,
\end{aligned}
```
with the update rule
```math
\begin{aligned}
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b^1_{i} V_{n,i}
                 + h \sum \limits_{i=1}^{s} b^2_{i} \Lambda_{n,i} , \\
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}^1_{i} F_{n,i}
                 + h \sum \limits_{i=1}^{s} \bar{b}^2_{i} G_{n,i} ,
\end{aligned}
```
and the definitions
```math
\begin{aligned}
F_{n,i} &= \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i})
         = - \nabla H (Q_{n,i}) + \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} , \\
G_{n,i} &= \nabla \vartheta (Q_{n,i}) \cdot \Lambda_{n,i} , \\
\Phi_{n,i} &= \phi (Q_{n,i}, P_{n,i}) = P_{n,i} - \vartheta (Q_{n,i}) , \\
\Psi_{n,i} &= \psi (Q_{n,i}, V_{n,i}, P_{n,i}, F_{n,i})
            = F_{n,i} - V_{n,i} \cdot \nabla \vartheta (Q_{n,i}) .
\end{aligned}
```

The four coefficient tableaus are stored as `q` ``\to a^1, b^1``, `q̃` ``\to a^2, b^2``,
`p` ``\to \bar{a}^1, \bar{b}^1`` and `p̃` ``\to \bar{a}^2, \bar{b}^2``. All
`SLRKLobattoIII*` constructors set `q̃ = q` and `p̃ = p`, so a method is fixed by one
Lobatto pair ``(a, \bar{a})``.

## The ``\omega`` constraint

``\omega`` is the ``s \times (s+1)`` matrix returned by `lobatto_ω_matrix(s)`. The
system ``\omega \, [\Psi_{n,1}, \ldots, \Psi_{n,s}, \phi(q_{n+1},p_{n+1})]^T = 0``
is equivalent to imposing the ``s-1`` Lobatto-IIIA-averaged secondary constraints
``\sum_j a^{\mathrm{IIIA}}_{ij} \Psi_{n,j} = 0`` for ``i = 2, ..., s`` **together
with** ``\phi (q_{n+1}, p_{n+1}) = 0``.

## The null vector

The Lobatto stage system is rank deficient by one in the ``V``-direction. The
multiplier ``\mu`` relaxes the primary constraint along the null vector ``d`` of the
Lagrange-derivative matrix (`get_lobatto_nullvector`), and the deficiency is removed
by the extra condition ``\sum_i d_i V_{n,i} = 0``. ``\mu`` enters the primary-constraint
equation and *only* that one, as in `VPARK`/`VSPARK`/`SPARK`. (`VPRK` implements the same
idea with a different residual layout: it has no primary-constraint row, so it perturbs
the momentum-stage equation instead — see `residual_correction!`.)

# Requirements on the problem

!!! warning
    `SLRK` requires the [`LDAEProblem`](@ref)'s `f` to be the **full** force
    ``\partial L / \partial q = - \nabla H + \nabla \vartheta \cdot v``, because `F`
    is used directly both in the ``p`` update and as ``\dot{p}`` in the secondary
    constraint ``\psi``. This is the *opposite* convention to `VSPARKsecondary`,
    which expects `f = -∇H` and reconstructs ``\partial L / \partial q`` as
    `f + g(V)`. Passing an `LDAEProblem` built for `VSPARKsecondary` to `SLRK`
    silently drops the ``\nabla \vartheta \cdot v`` term.

    `GeometricProblems.LotkaVolterra2d` ships the two variants side by side:
    use `ldaeproblem_slrk` with `SLRK`, and `ldaeproblem` with `VSPARKsecondary`.

# Properties

* The primary constraint ``\phi(q_n, p_n) = 0`` is preserved to round-off at every
  step (measured: ``\max |\phi| \approx 2 \times 10^{-15}`` over ``10^3`` steps on
  Lotka–Volterra), which is the point of the method.
* The order is ``2s-2`` for all six `SLRKLobattoIII*` families, confirmed empirically at
  ``s = 2, 3`` in `test/verification/spark_convergence_tests.jl`.
* **Symplecticity is approximate, not exact.** The underlying manuscript claims exact
  preservation of the noncanonical symplectic form, but its proof leaves an
  uncontrolled term in the multiplier block: constraint ``\sum_i d_i V_{n,i} = 0``
  has no counterpart for ``\Lambda``, so a residual
  ``h \sum_j b^2_j (d_j / \bar{b}^1_j) \, \mathrm{d}\mu \wedge \mathrm{d}\Lambda_{n,j}``
  survives. Measured on Lotka–Volterra, the first Poincaré invariant ``\oint p \, dq``
  drifts secularly at ``O(h^{p+1})`` per step, whereas a genuinely symplectic variational
  integrator holds the same (canonical) invariant to round-off out to ``3 \times 10^3``
  steps. See the "Fifth pass" section of the verification report,
  `docs/src/audit.md`, for the numbers.
* The defect is **gauge invariant**: two Lagrangians whose one-forms differ by an exact
  form give the same trajectory and the same drift, to round-off.

!!! warning "Choice of gauge"
    `SLRKLobattoIIIAB` and `SLRKLobattoIIIBA` are the two families for which **both**
    ``a^1`` and ``\bar{a}^1`` are rank deficient (Lobatto IIIA has a zero first row,
    IIIB a zero last column, and the conjugate partners inherit it). If a component of
    ``\vartheta`` vanishes identically, ``\nabla \vartheta^{T}`` loses rank, the
    projection force ``G = \nabla\vartheta^{T}\Lambda`` can no longer reach the
    corresponding momentum component, and their stage system is singular at *every*
    step size — e.g. they fail on
    `GeometricProblems.LotkaVolterra2dSingular` while running fine on the
    gauge-equivalent `LotkaVolterra2d`. `SLRKLobattoIIID` and `SLRKLobattoIIIE` are
    full rank in both blocks and are the safest default.
"""
const IntegratorSLRK{DT,TT} = GeometricIntegrator{<:LDAEProblem{DT,TT},<:SLRK}


function Base.show(io::IO, int::IntegratorSLRK)
    print(io, "\nSpecialised Partitioned Additive Runge-Kutta integrator for degenerate")
    print(io, "\nvariational systems with projection on secondary constraint:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(method(int)))\n")
    print(io, "   $(string(method(int).q))")
    print(io, "   $(string(method(int).p))")
    # print(io, reference(method(int)))
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:SLRK,<:LDAEProblem})
    # get caches for nonlinear solver vector
    local x = cache(int).x
    local D = ndims(cache(int))

    # compute initial guess for internal stages
    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        # Use the same node as `components!` (`q.c`); all Lobatto pairs share nodes.
        soltmp = (
            t=history[1].t + timestep(int) * tableau(int).q.c[i],
            q=cache(int).Qp[i],
            p=cache(int).Pp[i],
            q̇=cache(int).Vp[i],
            ṗ=cache(int).Fp[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in 1:pstages(method(int))
        for k in 1:D
            offset = 4 * (D * (i - 1) + k - 1)
            x[offset+1] = (cache(int).Qp[i][k] - sol.q[k]) / timestep(int)
            x[offset+2] = (cache(int).Pp[i][k] - sol.p[k]) / timestep(int)
            x[offset+3] = cache(int).Vp[i][k]
            x[offset+4] = 0
        end
    end

    if hasnullvector(method(int))
        for k in 1:D
            x[4*D*pstages(method(int))+k] = 0
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:SLRK}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)
    local D = ndims(C)

    for i in eachindex(C.Yp, C.Zp, C.Vp, C.Λp)
        for k in eachindex(C.Yp[i], C.Zp[i], C.Vp[i], C.Λp[i])
            # copy y to Y, Z and Λ
            C.Yp[i][k] = x[4*(D*(i-1)+k-1)+1]
            C.Zp[i][k] = x[4*(D*(i-1)+k-1)+2]
            C.Vp[i][k] = x[4*(D*(i-1)+k-1)+3]
            C.Λp[i][k] = x[4*(D*(i-1)+k-1)+4]

            # compute Q and P
            C.Qp[i][k] = sol.q[k] + timestep(int) * C.Yp[i][k]
            C.Pp[i][k] = sol.p[k] + timestep(int) * C.Zp[i][k]
        end

        # compute f(X)
        # f, g, ϕ and ψ are all evaluated on the position stage Qᵢ, which is generated
        # by the `q` tableau, so its node is `q.c[i]`. For every Lobatto tableau the
        # nodes of the pair coincide, so this agrees with `p.c[i]` for all shipped
        # `SLRKLobattoIII*` methods; `initial_guess!` uses the same tableau.
        tᵢ = sol.t + timestep(int) * (method(int).q.c[i] - 1)
        equations(int).f(C.Fp[i], tᵢ, C.Qp[i], C.Vp[i], params)
        equations(int).g(C.Gp[i], tᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Λp[i], params)
        equations(int).ϕ(C.Φp[i], tᵢ, C.Qp[i], C.Vp[i], C.Pp[i], params)
        equations(int).ψ(C.Ψp[i], tᵢ, C.Qp[i], C.Vp[i], C.Pp[i], C.Vp[i], C.Fp[i], params)
    end

    if hasnullvector(method(int))
        for k in eachindex(C.μ)
            C.μ[k] = x[4*D*nstages(int)+k]
        end
    end

    # compute q and p
    C.q̃ .= sol.q
    C.p̃ .= sol.p
    for i in 1:nstages(int)
        C.q̃ .+= timestep(int) .* method(int).q.b[i] .* C.Vp[i]
        C.q̃ .+= timestep(int) .* method(int).q̃.b[i] .* C.Λp[i]
        C.p̃ .+= timestep(int) .* method(int).p.b[i] .* C.Fp[i]
        C.p̃ .+= timestep(int) .* method(int).p̃.b[i] .* C.Gp[i]
    end

    # compute ϕ(q,p) at the new solution, for the last column of the ω constraint.
    # No velocity is available at (q_{n+1}, p_{n+1}); ϕ of a degenerate Lagrangian,
    # ϕ(q,p) = p - ϑ(q), does not use it. Zero it explicitly rather than relying on
    # the cache never having been written.
    C.ṽ .= 0
    equations(int).ϕ(C.ϕ̃, sol.t, C.q̃, C.ṽ, C.p̃, params)
end


# Compute stages of specialised partitioned additive Runge-Kutta methods for variational systems.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:SLRK}) where {ST}
    # get cache and number of internal stages
    local C = cache(int, ST)
    local S = nstages(int)
    local D = ndims(C)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute b = - [(Y-AV-AU), (Z-AF-AG), Φ, ωΨ]
    for i in 1:nstages(int)
        for k in 1:D
            b[4*(D*(i-1)+k-1)+1] = -C.Yp[i][k]
            b[4*(D*(i-1)+k-1)+2] = -C.Zp[i][k]
            b[4*(D*(i-1)+k-1)+3] = -C.Φp[i][k]
            b[4*(D*(i-1)+k-1)+4] = method(int).ω[i, S+1] * C.ϕ̃[k]

            for j in 1:nstages(int)
                b[4*(D*(i-1)+k-1)+1] += method(int).q.a[i, j] * C.Vp[j][k]
                b[4*(D*(i-1)+k-1)+1] += method(int).q̃.a[i, j] * C.Λp[j][k]
                b[4*(D*(i-1)+k-1)+2] += method(int).p.a[i, j] * C.Fp[j][k]
                b[4*(D*(i-1)+k-1)+2] += method(int).p̃.a[i, j] * C.Gp[j][k]
                b[4*(D*(i-1)+k-1)+4] += method(int).ω[i, j] * C.Ψp[j][k]
            end
        end
    end

    if hasnullvector(method(int))
        # The Lobatto stage system is rank-deficient by one in the V-direction.
        # The multiplier μ relaxes the primary constraint along the null vector d,
        #     ϕ(Q_{n,i}, P_{n,i}) = (d_i / b̄_i) μ ,
        # and the rank deficiency is removed by the extra condition Σ_i d_i V_{n,i} = 0
        # below. `vpark`, `vspark` and `spark` use the same convention: they add μ to
        # their own Φ row and nothing else. (`VPRK` implements the same idea with a
        # different residual layout — it has no Φ row at all, so `residual_correction!`
        # in `src/integrators/vi/vprk_integrator.jl` perturbs the P-stage equation
        # instead. That row is in P-space, not Z-space, so it carries no factor h.)
        #
        # μ must appear in the constraint row only, NOT additionally in the
        # momentum-stage row: that row lives in Z-space (P = p + h·Z), so the same
        # coefficient there carries an extra factor h and the two contributions
        # combine to (1-h)·μ·d_i/b̄_i — which makes the Jacobian exactly singular
        # at Δt = 1 and ill-conditioned near it.
        for i in 1:nstages(int)
            for k in 1:D
                b[4*(D*(i-1)+k-1)+3] += C.μ[k] * method(int).d[i] / method(int).p.b[i]
            end
        end

        for k in 1:D
            b[4*D*S+k] = 0
            for i in 1:nstages(int)
                b[4*D*S+k] -= C.Vp[i][k] * method(int).d[i]
            end
        end
    end
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:SLRK,<:LDAEProblem}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol.q, cache(int, DT).Vp, method(int).q, timestep(int))
    update!(sol.q, cache(int, DT).Λp, method(int).q̃, timestep(int))
    update!(sol.p, cache(int, DT).Fp, method(int).p, timestep(int))
    update!(sol.p, cache(int, DT).Gp, method(int).p̃, timestep(int))
end

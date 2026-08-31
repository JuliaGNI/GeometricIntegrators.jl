@doc raw"""
Degenerate Variational Runge-Kutta (DVRK) method for degenerate Lagrangian
systems of the form
```math
L (q, \dot{q}) = \vartheta (q) \cdot \dot{q} - H (q) ,
\qquad
\vartheta_{\mu} = 0
\quad \text{for} \quad
\mu = d/2+1, \, ..., \, d ,
```
that is, Lagrangians linear in the velocities whose symplectic potential has
``d/2`` identically vanishing components. The internal stages read
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) = \vartheta (Q_{n,i}) , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, F_{n,j} , &
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) ,
\end{aligned}
```
and the update rules are
```math
\begin{aligned}
q^{\mu}_{n+1} &= q^{\mu}_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V^{\mu}_{n,i} , &
p^{\mu}_{n+1} &= p^{\mu}_{n} + h \sum \limits_{i=1}^{s} b_{i} \, F^{\mu}_{n,i} ,
&& \mu = 1, \, ..., \, d/2 , \\
&&
p^{\mu}_{n+1} &= \vartheta^{\mu} (q_{n+1}) , && \mu = 1, \, ..., \, d .
\end{aligned}
```
Note that the quadrature updates are imposed only on the first ``d/2`` components.
The remaining components of ``q_{n+1}`` carry no quadrature update and are
determined implicitly by the constraint ``p_{n+1} = \vartheta (q_{n+1})``, which is
also what defines ``p_{n+1}`` on output.

The method preserves the noncanonical symplectic two-form
``\omega = \mathrm{d} \vartheta`` provided that

  * the symplectic potential is given in a gauge in which ``\vartheta_{\mu} = 0``
    for ``\mu > d/2``, as assumed above,
  * the ``d/2 \times d/2`` matrix ``\partial \vartheta_{\mu} / \partial q^{\nu}``,
    for ``\mu \le d/2 < \nu``, is invertible,
  * the coefficient matrix ``A = (a_{ij})`` is invertible,
  * the coefficients satisfy the symplecticity condition
    ``b_{i} a_{ij} + b_{j} a_{ji} = b_{i} b_{j}``, and
  * the momentum is initialised consistently with ``p_{0} = \vartheta (q_{0})``.

The first two conditions are independent: the second says that the constraint
``p_{n+1} = \vartheta (q_{n+1})`` can be solved for the components of ``q_{n+1}``
that carry no quadrature update, and it is also what makes ``\omega``
nondegenerate, but it does not by itself imply that the remaining components of
``\vartheta`` vanish.

The Gauss-Legendre tableaus satisfy the two conditions on the coefficients for any
number of stages, and `DVRK(Gauss(s))` attains the full order ``2s`` on Lagrangians
of the above form. Applied to a Lagrangian outside this class — for instance the same
system written in a gauge in which no component of ``\vartheta`` vanishes, which
leaves the dynamics unchanged — the method remains convergent but the order drops
to ``s``.

A Degenerate Variational Runge-Kutta method is instantiated by either
passing a Runge-Kutta tableau or a Runge-Kutta method:
```
DVRK(tableau::Tableau; check_conditions = true)
DVRK(method::RKMethod; check_conditions = true)
```
The constructor checks the two tableau conditions and warns if either is
violated; pass `check_conditions = false` to suppress the warnings.
"""
struct DVRK{TT} <: DVIMethod
    tableau::TT

    function DVRK(tableau::TT; check_conditions = true) where {TT <: Tableau}
        if check_conditions
            if !RungeKutta.issymplectic(tableau)
                @warn "The tableau $(tableau.name) does not satisfy the symplecticity " *
                      "condition b_i a_ij + b_j a_ji = b_i b_j. The resulting DVRK " *
                      "method will not preserve the symplectic structure."
            end
            # Test invertibility by rank, not by `det`: the determinant is not
            # scale-invariant and shrinks rapidly with the number of stages, so
            # `det(a) ≈ 0` would flag perfectly well-conditioned tableaus
            # (det(Gauss(10).a) ≈ 1.5e-12 at cond(a) ≈ 1.2e2).
            # The rank is computed in Float64 so that extended-precision tableaus
            # (e.g. Gauss(BigFloat, s)) do not require a generic SVD implementation.
            if rank(Matrix{Float64}(tableau.a)) < tableau.s
                @warn "The coefficient matrix of the tableau $(tableau.name) is " *
                      "singular. The DVRK method requires an invertible coefficient " *
                      "matrix and may fail to be well defined."
            end
        end
        new{TT}(tableau)
    end
end

DVRK(method::RKMethod, args...; kwargs...) = DVRK(tableau(method); kwargs...)

GeometricBase.tableau(method::DVRK) = method.tableau
GeometricBase.order(method::DVRK) = order(tableau(method))
eachstage(method::DVRK) = eachstage(tableau(method))
isexplicit(method::DVRK) = false
isimplicit(method::DVRK) = true
issymmetric(method::DVRK) = RungeKutta.issymmetric(tableau(method))
issymplectic(method::DVRK) = RungeKutta.issymplectic(tableau(method))

# Unlike the other DVI methods, symplecticity of DVRK depends on the tableau, so
# the unconditional `issymplectic(::Type{<:DVIMethod}) = true` of dvi_common.jl
# cannot be answered for the bare type.
issymplectic(::Type{DVRK}) = missing

@doc raw"""
Degenerate Variational Runge-Kutta integrator cache.

### Fields

* `q`: initial guess of solution
* `v`: initial guess of vector field
* `θ`: initial guess of symplectic potential
* `f`: initial guess of force field
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Θ`: implicit function of internal stages
* `F`: vector field of implicit function
"""
struct DVRKCache{DT, S} <: IODEIntegratorCache{DT}
    x::Vector{DT}

    q̄::Vector{DT}
    p̄::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Θ::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function DVRKCache{DT, S}(ics) where {DT, S}
        D = length(vec(ics.q))
        check_dvi_dimension(D)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        new(zeros(DT, D * (S + 1)), zeros(DT, D), zeros(DT, D), zeros(DT, D),
            zeros(DT, D), zeros(DT, D), zeros(DT, D), Q, V, Θ, F)
    end
end

function reset!(cache::DVRKCache, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

nlsolution(cache::DVRKCache) = cache.x

"""
    check_dvrk_initial_conditions(problem)

Warn if the initial momentum of `problem` is not consistent with the primary
constraint `p₀ = ϑ(t₀, q₀, v₀)`. Symplecticity of [`DVRK`](@ref) is conditional on
this consistency, and an inconsistent `p₀` silently degrades the method.
"""
function check_dvrk_initial_conditions(problem::Union{IODEProblem, LODEProblem})
    ics = initial_conditions(problem)
    q₀, p₀ = vec(ics.q), vec(ics.p)
    v₀ = haskey(ics, :v) ? vec(ics.v) : zero(q₀)
    ϑ₀ = zero(p₀)
    functions(problem).ϑ(ϑ₀, initial_conditions(problem).t, q₀, v₀, parameters(problem))
    tol = sqrt(eps(eltype(p₀))) * max(one(eltype(p₀)), maximum(abs, ϑ₀))
    if maximum(abs, p₀ .- ϑ₀) > tol
        @warn "The initial momentum is not consistent with the primary constraint " *
              "p₀ = ϑ(q₀) (max deviation $(maximum(abs, p₀ .- ϑ₀))). DVRK is only " *
              "symplectic for consistent initial conditions."
    end
    return nothing
end

function Cache{ST}(problem::Union{IODEProblem, LODEProblem}, method::DVRK; kwargs...) where {ST}
    check_dvrk_initial_conditions(problem)
    DVRKCache{ST, nstages(tableau(method))}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::Union{IODEProblem, LODEProblem}, method::DVRK) = DVRKCache{
    ST, nstages(tableau(method))}

function Base.show(io::IO, int::GeometricIntegrator{<:DVRK})
    print(io, "\nRunge-Kutta Integrator for Degenerate Lagrangians with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(tableau(int)))
end

function initial_guess!(sol, history, params, int::GeometricIntegrator{<:DVRK})
    # set some local variables for convenience
    local D = length(cache(int).q)
    local x = nlsolution(int)

    # compute initial guess for internal stages
    for i in eachstage(int)
        soltmp = (
            t = sol.t + timestep(int) * (tableau(int).c[i] - 1),
            q = cache(int).Q[i],
            p = cache(int).Θ[i],
            q̇ = cache(int).V[i],
            ṗ = cache(int).F[i]
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end
    for i in eachindex(cache(int).V)
        for k in eachindex(cache(int).V[i])
            x[D * (i - 1) + k] = cache(int).V[i][k]
        end
    end

    # compute initial guess for solution
    soltmp = (
        t = sol.t,
        q = cache(int).q,
        p = cache(int).θ,
        q̇ = cache(int).v,
        ṗ = cache(int).f
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    for k in 1:D
        x[D * nstages(int) + k] = cache(int).q[k]
    end
end

function components!(x::Vector{ST}, sol, params, int::GeometricIntegrator{<:DVRK}) where {ST}
    # set some local variables for convenience and clarity
    local D = length(cache(int, ST).q)
    local S = nstages(tableau(int))

    # copy x to V
    for i in eachindex(cache(int, ST).V)
        for k in eachindex(cache(int, ST).V[i])
            cache(int, ST).V[i][k] = x[D * (i - 1) + k]
        end
    end

    # copy x to q
    for k in eachindex(cache(int, ST).q)
        cache(int, ST).q[k] = x[D * S + k]
    end

    # compute Q = q + Δt A V, Θ = ϑ(Q), F = f(Q,V)
    for i in eachindex(cache(int, ST).Q, cache(int, ST).F, cache(int, ST).Θ)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        for k in eachindex(cache(int, ST).Q[i])
            y1 = y2 = zero(ST)
            for j in eachindex(cache(int, ST).V)
                y1 += tableau(int).a[i, j] * cache(int, ST).V[j][k]
                y2 += tableau(int).â[i, j] * cache(int, ST).V[j][k]
            end
            cache(int, ST).Q[i][k] = sol.q[k] + timestep(int) * (y1 + y2)
        end
        equations(int).ϑ(
            cache(int, ST).Θ[i], tᵢ, cache(int, ST).Q[i], cache(int, ST).V[i], params)
        equations(int).f(
            cache(int, ST).F[i], tᵢ, cache(int, ST).Q[i], cache(int, ST).V[i], params)
    end

    # compute θ = ϑ(q_{n+1}), the momentum at the end of the step
    #
    # NOTE: `cache.v` is not part of the nonlinear solution; it holds the value
    # extrapolated in `initial_guess!` and stays fixed throughout the Newton
    # iteration. This is exact for the degenerate Lagrangians DVRK targets, whose
    # symplectic potential ϑ = ϑ(q) does not depend on the velocity, but it would
    # be inconsistent for a velocity-dependent ϑ. DVRK is not applicable to such
    # systems in any case.
    equations(int).ϑ(cache(int, ST).θ, sol.t, cache(int, ST).q, cache(int, ST).v, params)
end

# Compute stages of fully implicit Runge-Kutta methods.
function residual!(b::Vector{ST}, sol, params, int::GeometricIntegrator{<:DVRK}) where {ST}
    # set some local variables for convenience and clarity
    local D = length(cache(int, ST).q)
    local S = nstages(tableau(int))

    # compute b
    for i in eachindex(cache(int, ST).Θ)
        for k in eachindex(cache(int, ST).Θ[i])
            y1 = y2 = zero(ST)
            for j in eachindex(cache(int, ST).F)
                y1 += tableau(int).a[i, j] * cache(int, ST).F[j][k]
                y2 += tableau(int).â[i, j] * cache(int, ST).F[j][k]
            end
            b[D * (i - 1) + k] = cache(int, ST).Θ[i][k] - sol.p[k] -
                                 timestep(int) * (y1 + y2)
        end
    end
    for k in 1:div(D, 2)
        y1 = y2 = zero(ST)
        for j in eachindex(cache(int, ST).F)
            y1 += tableau(int).b[j] * cache(int, ST).F[j][k]
            y2 += tableau(int).b̂[j] * cache(int, ST).F[j][k]
        end
        b[D * S + k] = cache(int, ST).θ[k] - sol.p[k] - timestep(int) * (y1 + y2)
    end
    for k in 1:div(D, 2)
        y1 = y2 = zero(ST)
        for j in eachindex(cache(int, ST).V)
            y1 += tableau(int).b[j] * cache(int, ST).V[j][k]
            y2 += tableau(int).b̂[j] * cache(int, ST).V[j][k]
        end
        b[D * S + div(D, 2) + k] = cache(int, ST).q[k] - sol.q[k] -
                                   timestep(int) * (y1 + y2)
    end
end

function update!(sol, params, int::GeometricIntegrator{<:DVRK}, DT)
    # compute final update
    sol.q .= cache(int, DT).q
    sol.p .= cache(int, DT).θ
end

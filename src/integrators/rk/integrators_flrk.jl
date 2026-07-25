@doc raw"""
Formal Lagrangian Runge-Kutta Method.

Applies a Runge-Kutta method to the *formal Lagrangian* of a noncanonical
Hamiltonian system. Given an [`LODE`](@ref) or [`IODE`](@ref) whose vector field
``v`` is determined by
```math
\dot{q} = v(q) = \bar{\Omega}^{-T} (q) \, \nabla H (q) ,
\qquad
\bar{\Omega}_{ij} (q) = \vartheta_{i,j} (q) - \vartheta_{j,i} (q) ,
```
the formal Lagrangian introduces an adjoint variable ``p`` and reads
```math
\bar{L} (q, \dot{q}, p, \dot{p}) = p^{T} \dot{q} - p^{T} v(q) + \vartheta^{T} (q) \, v(q) - H(q) ,
```
where the last two (momentum-modified) terms are chosen so that ``p = \vartheta(q)``
solves the adjoint equations. Variation with respect to ``p`` returns the physical
equations of motion, and variation with respect to ``q`` gives the adjoint equations
```math
\dot{p} = \dfrac{\partial \bar{L}}{\partial q} = f(q, v) + \left( \dfrac{\partial v}{\partial q} \right)^{T} \big( \vartheta(q) - p \big) ,
\qquad
f(q, v) = \big( \nabla \vartheta (q) \big)^{T} v - \nabla H (q) .
```

### Discretisation

The position is advanced by an ordinary implicit Runge-Kutta method applied to
``\dot{q} = v(q)``,
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} ,
\end{aligned}
```
and the adjoint variable by the corresponding Runge-Kutta method applied to the
adjoint equations with ``\bar{a}_{ij} = a_{ij}``, ``\bar{b}_{i} = b_{i}``,
```math
\begin{aligned}
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, \dot{P}_{n,j} , &
\dot{P}_{n,i} &= F_{n,i} + J_{n,i}^{T} \big( \Theta_{n,i} - P_{n,i} \big) , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} b_{i} \, \dot{P}_{n,i} ,
\end{aligned}
```
where ``\Theta_{n,i} = \vartheta (Q_{n,i})``, ``F_{n,i} = f (Q_{n,i}, V_{n,i})`` and
``J_{n,i} = \partial v / \partial q (Q_{n,i})``. Because the adjoint equations are
*linear* in ``P``, these stage equations are not solved by Newton iteration but by
one linear solve of size ``d s \times d s``,
```math
\big( I + h \, (a \otimes J^{T}) \big) \, P = \mathbb{1} \otimes p_{n} + h \, (a \otimes I) \, \big( F + J^{T} \Theta \big) ,
```
performed after the position stages have converged. The manuscript notes that the
adjoint variables need not be solved for in practice, which is why this is a
post-hoc computation rather than part of the nonlinear residual.

### Requirements and caveats

* The tableau must be **Gauss-Legendre** (or at least have an invertible ``a``); the
  elimination of the discrete adjoint momenta relies on it.
* The problem's ``\bar{v}`` must depend on ``q`` only. For example
  `GeometricProblems.HarmonicOscillator.lodeproblem` has ``\bar{v}(t,q,p) = p_1``,
  for which ``\partial \bar{v} / \partial q`` vanishes identically and the method is
  meaningless; use `degenerate_lodeproblem` or
  `GeometricProblems.LotkaVolterra2d.lodeproblem` instead.
* **Symplecticity is not established.** The reference states explicitly that
  variational Runge-Kutta methods are *not* symplectic with respect to the
  noncanonical symplectic form, promises "pseudo-symplecticity" and then proves
  nothing; conjugate symplecticity is recorded as an open conjecture in one draft
  and as an unproven theorem in another. There is likewise no order theorem — the
  order is inherited from the underlying tableau and is verified numerically in
  `test/verification/flrk_convergence_tests.jl`. Accordingly
  `issymplectic` returns `missing`, and the energy is expected to drift
  slowly rather than be conserved.

### Reference

Michael Kraus.
Variational Integrators for Noncanonical Hamiltonian Systems.
Working manuscript.

### Constructors

```
FLRK(tableau::Tableau)
FLRK(method::RKMethod)
```
"""
struct FLRK{TT<:Tableau} <: LODEMethod
    tableau::TT

    function FLRK(tableau::TT) where {TT<:Tableau}
        new{TT}(tableau)
    end
end

FLRK(method::RKMethod, args...; kwargs...) = FLRK(tableau(method))

GeometricBase.tableau(method::FLRK, args...; kwargs...) = method.tableau
GeometricBase.order(method::FLRK) = RungeKutta.order(tableau(method))

@inline nstages(method::FLRK) = RungeKutta.nstages(tableau(method))
@inline eachstage(method::FLRK) = RungeKutta.eachstage(tableau(method))

isexplicit(::Union{FLRK,Type{<:FLRK}}) = false
isimplicit(::Union{FLRK,Type{<:FLRK}}) = true
issymmetric(method::FLRK) = RungeKutta.issymmetric(tableau(method))

# Deliberately `missing`, not `true`: see the caveats in the docstring.
issymplectic(::Union{FLRK,Type{<:FLRK}}) = missing
isenergypreserving(::Union{FLRK,Type{<:FLRK}}) = false

isiodemethod(::Union{FLRK,Type{<:FLRK}}) = true

default_solver(::FLRK) = Newton()
default_iguess(::FLRK) = HermiteExtrapolation()

solversize(problem::AbstractProblemIODE, method::FLRK) =
    length(vec(initial_conditions(problem).q)) * nstages(method)


function Base.show(io::IO, int::GeometricIntegrator{<:FLRK})
    print(io, "\nFormal Lagrangian Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
end


@doc raw"""
Formal Lagrangian Runge-Kutta integrator cache.

### Fields

* `x`: nonlinear solver solution vector
* `q̃`, `p̃`, `ṽ`, `f̃`: temporaries for the initial guess
* `Q`, `V`, `Y`: internal stages of the position, its vector field and `Y = h A V`
* `Θ`: one-form ``\vartheta`` at the internal stages
* `P`: adjoint momentum at the internal stages
* `F`: force ``f = (\nabla \vartheta)^{T} v - \nabla H`` at the internal stages
* `G`: ``J^{T} (\Theta - P)`` at the internal stages
* `Ṗ`: total adjoint force ``F + G``
* `J`: Jacobians ``\partial v / \partial q`` at the internal stages
* `A`, `δP`: system matrix and right-hand side of the linear adjoint solve
* `jac`: automatic-differentiation Jacobian of ``\bar{v}``
"""
struct FLRKCache{ST,D,S,JT} <: IODEIntegratorCache{ST}
    x::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}

    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}

    Θ::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    G::Vector{Vector{ST}}
    Ṗ::Vector{Vector{ST}}

    J::Vector{Matrix{ST}}
    A::Matrix{ST}
    δP::Vector{ST}

    jac::JT

    function FLRKCache{ST,D,S}(jac::JT) where {ST,D,S,JT}
        x = zeros(ST, D * S)

        q̃ = zeros(ST, D)
        p̃ = zeros(ST, D)
        ṽ = zeros(ST, D)
        f̃ = zeros(ST, D)

        Q = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)

        Θ = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        G = create_internal_stage_vector(ST, D, S)
        Ṗ = create_internal_stage_vector(ST, D, S)

        J = [zeros(ST, D, D) for _ in 1:S]
        A = zeros(ST, D * S, D * S)
        δP = zeros(ST, D * S)

        new{ST,D,S,JT}(x, q̃, p̃, ṽ, f̃, Q, V, Y, Θ, P, F, G, Ṗ, J, A, δP, jac)
    end
end

nlsolution(cache::FLRKCache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::FLRK; kwargs...) where {ST}
    D = length(vec(initial_conditions(problem).q))
    S = nstages(method)

    # ∂v̄/∂q at fixed (t, p, params). The extra arguments travel through the
    # `Jacobian` functor's `params` slot as a NamedTuple.
    local v̄ = initialguess(problem).v
    v_q!(v, q, p) = v̄(v, p.t, q, p.p, p.params)

    FLRKCache{ST,D,S}(Jacobian{ST}(v_q!, D, D); kwargs...)
end

# The Jacobian parameter is left free: the concrete `JacobianAutodiff` type depends on
# the problem's `v̄` closure, and type parameters are invariant, so naming it here
# would make the `CacheDict` type assertion fail.
@inline function CacheType(ST, problem::AbstractProblemIODE, method::FLRK)
    D = length(vec(initial_conditions(problem).q))
    FLRKCache{ST,D,nstages(method)}
end


function internal_variables(method::FLRK, problem::AbstractProblemIODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = length(vec(initial_conditions(problem).q))

    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    Θ = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)

    (Q=Q, V=V, Θ=Θ, P=P, F=F)
end

function copy_internal_variables!(solstep::SolutionStep, cache::FLRKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :Θ) && copyto!(internal(solstep).Θ, cache.Θ)
    haskey(internal(solstep), :P) && copyto!(internal(solstep).P, cache.P)
    haskey(internal(solstep), :F) && copyto!(internal(solstep).F, cache.F)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:FLRK})
    local x = nlsolution(int)
    local D = length(cache(int).q̃)

    # compute initial guess for the internal stages
    for i in eachstage(int)
        soltmp = (
            t=history[1].t + timestep(int) * tableau(int).c[i],
            q=cache(int).Q[i],
            p=cache(int).Θ[i],
            q̇=cache(int).V[i],
            ṗ=cache(int).F[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
        initialguess(problem(int)).v(soltmp.q̇, soltmp.t, soltmp.q, soltmp.p, params)
    end

    # assemble the solution vector, which holds Y = h A V
    for i in eachstage(int)
        offset = D * (i - 1)
        for k in 1:D
            x[offset+k] = 0
            for j in eachstage(int)
                x[offset+k] += timestep(int) * tableau(int).a[i, j] * cache(int).V[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:FLRK}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)

    # copy x to Y and compute Q = q + Y
    for i in eachindex(C.Q, C.Y)
        for k in eachindex(C.Q[i], C.Y[i])
            C.Y[i][k] = x[D*(i-1)+k]
            C.Q[i][k] = sol.q[k] + C.Y[i][k]
        end
    end

    # compute V = v̄(Q)
    for i in eachindex(C.Q, C.V)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        initialguess(problem(int)).v(C.V[i], tᵢ, C.Q[i], sol.p, params)
    end
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:FLRK}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)

    # compute b = Y - h A V
    for i in eachindex(C.Y)
        for k in eachindex(C.Y[i])
            y1 = y2 = zero(ST)
            for j in eachindex(C.V)
                y1 += tableau(int).a[i, j] * C.V[j][k]
                y2 += tableau(int).â[i, j] * C.V[j][k]
            end
            b[D*(i-1)+k] = C.Y[i][k] - timestep(int) * (y1 + y2)
        end
    end
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:FLRK}) where {ST}
    @assert axes(x) == axes(b)
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end


@doc raw"""
Solve the linear stage equations of the adjoint (momentum) variable and assemble the
total adjoint force ``\dot{P} = F + J^{T} (\Theta - P)``.

Runs on the `DT` cache only, after the position stages have converged.
"""
function adjoint_components!(sol, params, int::GeometricIntegrator{<:FLRK}, DT)
    local C = cache(int, DT)
    local D = length(C.q̃)
    local S = nstages(int)
    local h = timestep(int)

    # compute Θ = ϑ(Q,V), F = ∂L/∂q = f(Q,V) and J = ∂v̄/∂q (Q)
    for i in eachstage(int)
        tᵢ = sol.t + h * (tableau(int).c[i] - 1)
        equations(int).ϑ(C.Θ[i], tᵢ, C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], tᵢ, C.Q[i], C.V[i], params)
        C.jac(C.J[i], C.Q[i], (t=tᵢ, p=sol.p, params=params))
    end

    # G = Jᵀ Θ
    for l in eachstage(int)
        for i in 1:D
            C.G[l][i] = 0
            for j in 1:D
                C.G[l][i] += C.Θ[l][j] * C.J[l][j, i]
            end
        end
    end

    # δP = 1 ⊗ p + h (a ⊗ I) (F + Jᵀ Θ)
    for l in eachstage(int)
        for i in 1:D
            C.δP[(l-1)*D+i] = sol.p[i]
            for k in eachstage(int)
                C.δP[(l-1)*D+i] += h * tableau(int).a[l, k] * (C.F[k][i] + C.G[k][i])
            end
        end
    end

    # A = I + h (a ⊗ Jᵀ)
    fill!(C.A, 0)
    for k in eachstage(int)
        for l in eachstage(int)
            for i in 1:D
                for j in 1:D
                    C.A[(k-1)*D+i, (l-1)*D+j] = h * tableau(int).a[k, l] * C.J[l][j, i]
                end
            end
        end
    end
    for i in 1:(D*S)
        C.A[i, i] += one(DT)
    end

    # solve A P = δP
    LinearAlgebra.ldiv!(LinearAlgebra.lu!(C.A), C.δP)

    for l in eachstage(int)
        for i in 1:D
            C.P[l][i] = C.δP[(l-1)*D+i]
        end
    end

    # G = Jᵀ (Θ - P) and Ṗ = F + G
    for l in eachstage(int)
        for i in 1:D
            C.G[l][i] = 0
            for j in 1:D
                C.G[l][i] += (C.Θ[l][j] - C.P[l][j]) * C.J[l][j, i]
            end
            C.Ṗ[l][i] = C.F[l][i] + C.G[l][i]
        end
    end
end


function update!(sol, params, int::GeometricIntegrator{<:FLRK}, DT)
    # advance the position with the Runge-Kutta method applied to q̇ = v̄(q)
    update!(sol.q, cache(int, DT).V, tableau(int), timestep(int))

    # advance the adjoint momentum with the same tableau
    update!(sol.p, cache(int, DT).Ṗ, tableau(int), timestep(int))
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:FLRK}) where {DT}
    components!(x, sol, params, int)
    adjoint_components!(sol, params, int, DT)
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:FLRK,<:AbstractProblemIODE})
    # solve the position stages
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # compute the adjoint momentum and the final update
    update!(sol, params, nlsolution(int), int)
end

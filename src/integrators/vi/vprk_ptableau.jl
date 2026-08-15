@doc raw"""
Variational Partitioned Runge-Kutta method with projection **in the tableau**, for
degenerate Lagrangian systems. *EXPERIMENTAL.*

Where the other projected VPRK methods append a projection step to an unmodified
Runge-Kutta method, this one perturbs the tableau itself. It uses the family of
[`CoefficientsPGLRK`](@ref),
```math
a(\lambda) = a + A(\lambda) ,
\qquad
A(\lambda) = P \, W(\lambda) \, Q ,
```
where ``a`` is the ``s``-stage Gauss tableau and ``W(\lambda)`` is skew with the
``d`` free parameters ``\lambda`` placed in its high-index corner,
```math
W_{s-k+1,\,s-k} = + \lambda_{k} ,
\qquad
W_{s-k,\,s-k+1} = - \lambda_{k} ,
\qquad
k = 1, \dots, d .
```
The stage equations are those of a VPRK method with the perturbed tableau,
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} \big( a_{ij} + A_{ij} (\lambda) \big) \, V_{n,j} , &
P_{n,i} &= \vartheta (Q_{n,i}, V_{n,i}) , \\
Z_{n,i} &= \sum \limits_{j=1}^{s} \big( a_{ij} + A_{ij} (\lambda) \big) \, F_{n,j} , &
P_{n,i} &= p_{n} + h \, Z_{n,i} ,
\end{aligned}
```
and ``\lambda`` is fixed by requiring that the **Dirac constraint** hold at the end of
the step,
```math
\vartheta (q_{n+1}) - p_{n+1} = 0 ,
\qquad
q_{n+1} = q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} ,
\qquad
p_{n+1} = p_{n} + h \sum \limits_{i=1}^{s} b_{i} \, F_{n,i} .
```
That is ``d`` scalar conditions for the ``d`` components of ``\lambda``, so the stage
equations and the constraint are solved together as one nonlinear system of size
``d (s+1)``.

Because ``B A(\lambda)`` is skew for ``B = \mathrm{diag}(b)``, the perturbed tableau
satisfies the symplecticity condition for every ``\lambda``.

### Stage-count requirement

The ``k``-th multiplier occupies the skew slot ``(s-k+1, s-k)``, so ``d`` multipliers
need ``s \geq d + 1`` merely to fit into ``W`` — this is asserted. Full order
additionally requires ``s \geq d + 2``, because the slot ``(2,1)`` is the one pair that
coincides with the order-determining ``\xi_1`` entries of ``X`` and destroys the
conditions ``b^{T} A = 0`` and ``A \mathbb{1} = 0``. At ``s = d + 1`` the method
therefore remains symplectic but loses order. For a two-dimensional problem that means
``s \geq 4`` for full order.

### Status

This method has **no written-up reference**. The nearest documented relative,
`sec:glrk` of *Notes on Projected VIs for Degenerate Lagrangian Systems*, places the
multiplier in the stage equations rather than in the tableau, and its symplecticity
computation is truncated; the *Internal Projection* section of
*Projected Variational Integrators for Degenerate Lagrangian Systems* that would cover
this construction is an empty stub. The tableau family itself is that of
Brugnano, Iavernaro & Trigiante (see [`CoefficientsPGLRK`](@ref)), with the Dirac
constraint substituted for their energy constraint. Accordingly the order is measured
rather than claimed, and `issymplectic` reflects only the algebraic property
above.

### Constructors

```
VPRKpTableau(coefficients::CoefficientsPGLRK)
VPRKpTableau(s::Int, T = Float64)
```
"""
struct VPRKpTableau{CT<:CoefficientsPGLRK} <: VIMethod
    coefficients::CT
end

VPRKpTableau(s::Int, ::Type{T}=Float64) where {T} = VPRKpTableau(CoefficientsPGLRK(T, s))

GeometricBase.tableau(method::VPRKpTableau, args...; kwargs...) = method.coefficients
GeometricBase.order(method::VPRKpTableau) = RungeKutta.order(tableau(method))
GeometricBase.order(::Type{<:VPRKpTableau}) = "2s"

@inline nstages(method::VPRKpTableau) = RungeKutta.nstages(tableau(method))
@inline eachstage(method::VPRKpTableau) = RungeKutta.eachstage(tableau(method))

isexplicit(::Union{VPRKpTableau,Type{<:VPRKpTableau}}) = false
isimplicit(::Union{VPRKpTableau,Type{<:VPRKpTableau}}) = true
issymmetric(::Union{VPRKpTableau,Type{<:VPRKpTableau}}) = missing

# `B A(λ)` is skew for every λ, so the perturbed tableau satisfies the symplecticity
# condition. Whether the *projected* method is symplectic is not established: the
# relevant derivation in the reference is truncated.
issymplectic(::Union{VPRKpTableau,Type{<:VPRKpTableau}}) = missing

# stage equations (D*S) plus the Dirac constraint (D)
solversize(method::VPRKpTableau, problem::AbstractProblemIODE) =
    length(vec(initial_conditions(problem).q)) * (nstages(method) + 1)


function Base.show(io::IO, int::GeometricIntegrator{<:VPRKpTableau})
    print(io, "\nVariational Partitioned Runge-Kutta Integrator with Projection in the Tableau with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
end


@doc raw"""
Cache of [`VPRKpTableau`](@ref).

`W`, `T` and `A` are `ST`-typed because ``\lambda`` is one of the nonlinear solver
unknowns, so the perturbed tableau carries dual numbers during the Jacobian
computation.

### Fields

* `x`: nonlinear solver solution vector, holding `(V₁, …, V_S, λ)`
* `λ`: the Dirac-constraint multipliers
* `q`, `p`, `θ`: solution and one-form at the end of the step
* `z`: permanent zeros, for the unused velocity slot of `ϑ`
* `q̃`, `p̃`, `ṽ`, `f̃`: temporaries for the initial guess
* `Q`, `P`, `V`, `F`, `Y`, `Z`: internal stages
* `W`, `T`, `A`: skew generator, workspace and the perturbation ``A = P W Q``
"""
struct VPRKpTableauCache{ST,D,S} <: IODEIntegratorCache{ST}
    x::Vector{ST}

    λ::Vector{ST}

    q::Vector{ST}
    p::Vector{ST}
    θ::Vector{ST}

    # permanent zeros, for the unused velocity slot of `ϑ`; never written
    z::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}
    Z::Vector{Vector{ST}}

    W::Matrix{ST}
    T::Matrix{ST}
    A::Matrix{ST}

    function VPRKpTableauCache{ST,D,S}() where {ST,D,S}
        x = zeros(ST, D * (S + 1))

        λ = zeros(ST, D)

        q = zeros(ST, D)
        p = zeros(ST, D)
        θ = zeros(ST, D)
        z = zeros(ST, D)

        q̃ = zeros(ST, D)
        p̃ = zeros(ST, D)
        ṽ = zeros(ST, D)
        f̃ = zeros(ST, D)

        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        W = zeros(ST, S, S)
        T = zeros(ST, S, S)
        A = zeros(ST, S, S)

        new(x, λ, q, p, θ, z, q̃, p̃, ṽ, f̃, Q, P, V, F, Y, Z, W, T, A)
    end
end

nlsolution(cache::VPRKpTableauCache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::VPRKpTableau; kwargs...) where {ST}
    D = length(vec(initial_conditions(problem).q))
    S = nstages(method)
    @assert S ≥ D + 1 "VPRKpTableau needs at least D+1 = $(D+1) stages to fit $(D) " *
                      "multipliers into the skew generator W, and s ≥ D+2 = $(D+2) " *
                      "stages for full order (got s = $(S))."
    VPRKpTableauCache{ST,D,S}(; kwargs...)
end

@inline function CacheType(ST, problem::AbstractProblemIODE, method::VPRKpTableau)
    D = length(vec(initial_conditions(problem).q))
    VPRKpTableauCache{ST,D,nstages(method)}
end


function internal_variables(method::VPRKpTableau, problem::AbstractProblemIODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = length(vec(initial_conditions(problem).q))

    Q = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)

    (Q=Q, P=P, V=V, F=F)
end

function copy_internal_variables!(solstep::SolutionStep, cache::VPRKpTableauCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :P) && copyto!(internal(solstep).P, cache.P)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :F) && copyto!(internal(solstep).F, cache.F)
end


"""
Recompute the tableau perturbation `A = P W(λ) Q` for the current multipliers.

The `k`-th multiplier is placed in the skew slot `(s-k+1, s-k)`, i.e. as far into the
high-index corner of `W` as possible, which is where the perturbation leaves the order
conditions intact.
"""
function update_tableau!(cache::VPRKpTableauCache, coefficients::CoefficientsPGLRK)
    local S = size(cache.W, 1)

    fill!(cache.W, 0)
    for k in eachindex(cache.λ)
        cache.W[S-k+1, S-k] = +cache.λ[k]
        cache.W[S-k, S-k+1] = -cache.λ[k]
    end

    mul!(cache.T, cache.W, coefficients.Q)
    mul!(cache.A, coefficients.P, cache.T)

    return cache.A
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:VPRKpTableau})
    local x = nlsolution(int)
    local D = length(cache(int).q̃)
    local S = nstages(int)

    for i in eachstage(int)
        soltmp = (
            t=history[1].t + timestep(int) * tableau(int).c[i],
            q=cache(int).Q[i],
            p=cache(int).P[i],
            q̇=cache(int).V[i],
            ṗ=cache(int).F[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:D
            x[D*(i-1)+k] = cache(int).V[i][k]
        end
    end

    # start from the unperturbed tableau
    for k in 1:D
        x[D*S+k] = 0
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRKpTableau}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)
    local S = nstages(int)

    # copy x to V and λ
    for i in eachindex(C.V)
        for k in eachindex(C.V[i])
            C.V[i][k] = x[D*(i-1)+k]
        end
    end
    for k in eachindex(C.λ)
        C.λ[k] = x[D*S+k]
    end

    # rebuild the perturbed tableau for the current λ
    update_tableau!(C, tableau(int))

    # compute Y = (a + A) V and Q = q + Δt Y
    for i in eachindex(C.Q, C.Y)
        for k in eachindex(C.Q[i])
            y = zero(ST)
            for j in eachindex(C.V)
                y += (tableau(int).a[i, j] + C.A[i, j]) * C.V[j][k]
            end
            C.Y[i][k] = y
            C.Q[i][k] = sol.q[k] + timestep(int) * y
        end
    end

    # compute P = ϑ(Q,V) and F = f(Q,V)
    for i in eachindex(C.P, C.F)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        equations(int).ϑ(C.P[i], tᵢ, C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], tᵢ, C.Q[i], C.V[i], params)
    end

    # compute Z = (a + A) F
    for i in eachindex(C.Z)
        for k in eachindex(C.Z[i])
            z = zero(ST)
            for j in eachindex(C.F)
                z += (tableau(int).a[i, j] + C.A[i, j]) * C.F[j][k]
            end
            C.Z[i][k] = z
        end
    end

    # compute the solution at the end of the step
    for k in eachindex(C.q, C.p)
        y = z = zero(ST)
        for j in eachindex(C.V, C.F)
            y += tableau(int).b[j] * C.V[j][k]
            z += tableau(int).b[j] * C.F[j][k]
        end
        C.q[k] = sol.q[k] + timestep(int) * y
        C.p[k] = sol.p[k] + timestep(int) * z
    end

    # compute θ = ϑ(q) for the Dirac constraint. `ϑ` needs a velocity slot, which is
    # supplied as the permanent zeros `C.z`: this method targets degenerate Lagrangians,
    # for which `ϑ` does not depend on the velocity. On a `ϑ` that *does* read it the
    # constraint would silently be imposed at v = 0.
    equations(int).ϑ(C.θ, sol.t, C.q, C.z, params)
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRKpTableau}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)
    local S = nstages(int)

    # stage equations: P - p - Δt Z = 0
    for i in eachindex(C.P, C.Z)
        for k in eachindex(C.P[i])
            b[D*(i-1)+k] = C.P[i][k] - sol.p[k] - timestep(int) * C.Z[i][k]
        end
    end

    # Dirac constraint: ϑ(q_{n+1}) - p_{n+1} = 0
    for k in eachindex(C.θ, C.p)
        b[D*S+k] = C.θ[k] - C.p[k]
    end
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRKpTableau}) where {ST}
    @assert axes(x) == axes(b)
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end


function update!(sol, params, int::GeometricIntegrator{<:VPRKpTableau}, DT)
    sol.q .= cache(int, DT).q
    sol.p .= cache(int, DT).p
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:VPRKpTableau}) where {DT}
    components!(x, sol, params, int)
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:VPRKpTableau,<:AbstractProblemIODE})
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))
    check_solver_status(solverstatus, int)
    update!(sol, params, nlsolution(int), int)
end

@doc raw"""
Degenerate Variational Runge-Kutta (DVRK) method for noncanonical
symplectic equations solving the system
```math
\begin{aligned}
P_{n,i} &= \vartheta (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} ,
\end{aligned}
```
Usually we are interested in Lagrangian systems, where
```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) ,
\end{aligned}
```
and tableaus satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```

A Degenerate Variational Runge-Kutta method is instantiated by either
passing a Runge-Kutta tableau or a Runge-Kutta method:
```
DVRK(tableau::Tableau)
DVRK(method::RKMethod)
```
"""
struct DVRK{TT} <: DVIMethod
    tableau::TT

    function DVRK(tableau::TT) where {TT <: Tableau}
        new{TT}(tableau)
    end
end

DVRK(method::RKMethod, args...; kwargs...) = DVRK(tableau(method))

GeometricBase.tableau(method::DVRK) = method.tableau
GeometricBase.order(method::DVRK) = order(tableaus(method))
eachstage(method::DVRK) = eachstage(tableau(method))
isexplicit(method::DVRK) = false
isimplicit(method::DVRK) = true
issymmetric(method::DVRK) = issymmetric(tableaus(method))
issymplectic(method::DVRK) = issymplectic(tableaus(method))


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
struct DVRKCache{DT,D,S} <: IODEIntegratorCache{DT,D}
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

    function DVRKCache{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        new(zeros(DT, D*(S+1)), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Θ, F)
    end
end

function reset!(cache::DVRKCache, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

nlsolution(cache::DVRKCache) = cache.x

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::DVRK; kwargs...) where {ST}
    DVRKCache{ST, ndims(problem), nstages(tableau(method))}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::DVRK) = DVRKCache{ST, ndims(problem)}


function Base.show(io::IO, int::GeometricIntegrator{<:DVRK})
    print(io, "\nRunge-Kutta Integrator for Degenerate Lagrangians with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(tableau(int)))
end


function initial_guess!(int::GeometricIntegrator{<:DVRK})
    # set some local variables for convenience
    local D = ndims(int)
    local x = nlsolution(int)

    # compute initial guess for internal stages
    for i in eachstage(int)
        initialguess!(solstep(int).t̄ + timestep(int) * tableau(int).c[i], cache(int).Q[i], cache(int).Θ[i], cache(int).V[i], cache(int).F[i], solstep(int), problem(int), iguess(int))
    end
    for i in eachindex(cache(int).V)
        for k in eachindex(cache(int).V[i])
            x[D*(i-1)+k] = cache(int).V[i][k]
        end
    end
    
    # compute initial guess for solution
    initialguess!(solstep(int).t, cache(int).q, cache(int).θ, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

    for k in 1:D
        x[ndims(int) * nstages(int) + k] = cache(int).q[k]
    end
end


function components!(x::Vector{ST}, int::GeometricIntegrator{<:DVRK}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local S = nstages(tableau(int))

    # copy x to V
    for i in eachindex(cache(int,ST).V)
        for k in eachindex(cache(int,ST).V[i])
            cache(int,ST).V[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to q
    for k in eachindex(cache(int,ST).q)
        cache(int,ST).q[k] = x[D*S+k]
    end

    # compute Q = q + Δt A V, Θ = ϑ(Q), F = f(Q,V)
    for i in eachindex(cache(int,ST).Q, cache(int,ST).F, cache(int,ST).Θ)
        tᵢ = solstep(int).t̄ + timestep(int) * tableau(int).c[i]
        for k in eachindex(cache(int,ST).Q[i])
            y1 = y2 = zero(ST)
            for j in eachindex(cache(int,ST).V)
                y1 += tableau(int).a[i,j] * cache(int,ST).V[j][k]
                y2 += tableau(int).â[i,j] * cache(int,ST).V[j][k]
            end
            cache(int,ST).Q[i][k] = cache(int).q̄[k] + timestep(int) * (y1 + y2)
        end
        equations(int).ϑ(cache(int,ST).Θ[i], tᵢ, cache(int,ST).Q[i], cache(int,ST).V[i], parameters(solstep(int)))
        equations(int).f(cache(int,ST).F[i], tᵢ, cache(int,ST).Q[i], cache(int,ST).V[i], parameters(solstep(int)))
    end

    # compute q̄ = q + Δt B V, Θ = ϑ(q̄)
    equations(int).ϑ(cache(int,ST).θ, solstep(int).t, cache(int,ST).q, cache(int,ST).v, parameters(solstep(int)))
end


# Compute stages of fully implicit Runge-Kutta methods.
function residual!(b::Vector{ST}, x::Vector{ST}, int::GeometricIntegrator{<:DVRK}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local S = nstages(tableau(int))

    # compute stages from nonlinear solver solution x
    components!(x, int)
    
    # compute b
    for i in eachindex(cache(int,ST).Θ)
        for k in eachindex(cache(int,ST).Θ[i])
            y1 = y2 = zero(ST)
            for j in eachindex(cache(int,ST).F)
                y1 += tableau(int).a[i,j] * cache(int,ST).F[j][k]
                y2 += tableau(int).â[i,j] * cache(int,ST).F[j][k]
            end
            b[D*(i-1)+k] = cache(int,ST).Θ[i][k] - cache(int).p̄[k] - timestep(int) * (y1 + y2)
        end
    end
    for k in 1:div(D,2)
        y1 = y2 = zero(ST)
        for j in eachindex(cache(int,ST).F)
            y1 += tableau(int).b[j] * cache(int,ST).F[j][k]
            y2 += tableau(int).b̂[j] * cache(int,ST).F[j][k]
        end
        b[D*S+k] = cache(int,ST).θ[k] - cache(int).p̄[k] - timestep(int) * (y1 + y2)
    end
    for k in 1:div(D,2)
        y1 = y2 = zero(ST)
        for j in eachindex(cache(int,ST).V)
            y1 += tableau(int).b[j] * cache(int,ST).V[j][k]
            y2 += tableau(int).b̂[j] * cache(int,ST).V[j][k]
        end
        b[D*S+div(D,2)+k] = cache(int,ST).q[k] - cache(int).q̄[k] - timestep(int) * (y1 + y2)
    end
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:DVRK}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int, DT).q
    solstep(int).p .= cache(int, DT).θ
end

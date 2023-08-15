@doc raw"""
Implicit Runge-Kutta integrator cache.

### Fields

* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Θ`: internal stages of one-form ``\vartheta``
* `F`: internal stages of force field
"""
struct IntegratorCacheIRKimplicit{DT,D,S} <: IODEIntegratorCache{DT,D}
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

    function IntegratorCacheIRKimplicit{DT,D,S}(method::IRKMethod) where {DT,D,S}
        if implicit_update(method)
            x = zeros(DT, D*(S+1))
        else
            x = zeros(DT, D*S)
        end

        q̄ = zeros(DT,D)
        p̄ = zeros(DT,D)

        q = zeros(DT,D)
        v = zeros(DT,D)
        θ = zeros(DT,D)
        f = zeros(DT,D)

        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        
        new(x, q̄, p̄, q, v, θ, f, Q, V, Θ, F)
    end
end

nlsolution(cache::IntegratorCacheIRKimplicit) = cache.x

function reset!(cache::IntegratorCacheIRKimplicit, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::IRKMethod; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheIRKimplicit{ST,D,S}(method; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::IRKMethod) = IntegratorCacheIRKimplicit{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Implicit Runge-Kutta integrator for implicit systems of equations,
solving the system
```math
\begin{aligned}
P_{n,i} &= \vartheta (Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , \\
F_{n,i} &= f (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} .
\end{aligned}
```
If `implicit_update` is set to `true`, the update is computed by solving
```math
\vartheta(q_{n+1}) = \vartheta(q_{n}) + h \sum \limits_{i=1}^{s} b_{i}  \, f (Q_{n,j}, V_{n,j}) ,
```
otherwise it is computed explicitly by
```math
q_{n+1} = q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
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
"""
const IntegratorIRKimplicit{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:IRKMethod}

function solversize(problem::Union{IODEProblem,LODEProblem}, method::IRKMethod)
    n = ndims(problem) * nstages(method)

    if implicit_update(method)
        n += ndims(problem)
    end

    return n
end


function Base.show(io::IO, int::IntegratorIRKimplicit)
    print(io, "\nRunge-Kutta Integrator for Implicit Equations with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


function initial_guess!(
    solstep::SolutionStepPODE{DT},
    problem::Union{IODEProblem,LODEProblem},
    method::IRKMethod,
    caches::CacheDict,
    ::NonlinearSolver,
    iguess::Union{InitialGuess,Extrapolation}) where {DT}
    
    local cache = caches[DT]
    local offset::Int

    # compute initial guess for internal stages
    for i in eachstage(method)
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).c[i], cache.Q[i], cache.Θ[i], cache.V[i], cache.F[i], solstep, problem, iguess)
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(method)
        offset = ndims(problem)*(i-1)
        for k in 1:ndims(problem)
            cache.x[offset+k] = cache.Θ[i][k] - solstep.p̄[k]
            for j in eachstage(method)
                cache.x[offset+k] -= timestep(problem) * tableau(method).a[i,j] * cache.F[j][k]
            end
        end
    end

    # compute initial guess for solution
    if implicit_update(method)
        initialguess!(solstep.t, cache.q, cache.θ, cache.v, cache.f, solstep, problem, iguess)

        offset = ndims(problem) * nstages(tableau(method))

        for k in 1:ndims(problem)
            cache.x[offset + k] = cache.q[k]
        end
    end
end


function components!(x::AbstractVector{ST}, int::IntegratorIRKimplicit) where {ST}
    # temporary variables
    local y1::ST
    local y2::ST
    local tᵢ::timetype(problem(int))

    # get cache for internal stages
    local q̄ = cache(int, ST).q̄
    local q = cache(int, ST).q
    local v = cache(int, ST).v
    local θ = cache(int, ST).θ
    local Q = cache(int, ST).Q
    local V = cache(int, ST).V
    local Θ = cache(int, ST).Θ
    local F = cache(int, ST).F

    # copy x to V
    for i in eachindex(V)
        for k in eachindex(V[i])
            V[i][k] = x[ndims(int)*(i-1) + k]
        end
    end

    # copy x to q
    if implicit_update(int)
        for k in eachindex(q)
            q[k] = x[ndims(int) * nstages(int) + k]
        end
    end

    # compute Q = q + Δt A V
    for i in eachindex(Q)
        for k in eachindex(Q[i])
            y1 = y2 = 0
            for j in eachindex(V)
                y1 += tableau(int).a[i,j] * V[j][k]
                y2 += tableau(int).â[i,j] * V[j][k]
            end
            Q[i][k] = q̄[k] + timestep(int) * (y1 + y2)
        end
    end

    # compute Θ = ϑ(Q) and F = f(Q,V)
    for i in eachindex(Θ,F)
        tᵢ = solstep(int).t̄ + timestep(int) * tableau(int).c[i]
        equations(int).ϑ(Θ[i], tᵢ, Q[i], V[i])
        equations(int).f(F[i], tᵢ, Q[i], V[i])
    end

    # compute θ = ϑ(q)
    if implicit_update(int)
        equations(int).ϑ(θ, solstep(int).t, q, v)
    end
end


# Compute stages of implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, int::IntegratorIRKimplicit) where {ST}
    # get cache for previous solution and internal stages
    local p̄ = cache(int, ST).p̄
    local θ = cache(int, ST).θ
    local Θ = cache(int, ST).Θ
    local F = cache(int, ST).F

    # temporary variables
    local y1::ST
    local y2::ST

    # compute b for internal stages
    for i in eachindex(Θ)
        for k in eachindex(Θ[i])
            y1 = y2 = 0
            for j in eachindex(F)
                y1 += tableau(int).a[i,j] * F[j][k]
                y2 += tableau(int).â[i,j] * F[j][k]
            end
            b[ndims(int)*(i-1) + k] = Θ[i][k] - p̄[k] - timestep(int) * (y1 + y2)
        end
    end

    # compute b for update
    if implicit_update(int)
        for k in eachindex(θ)
            y1 = y2 = 0
            for j in eachindex(F)
                y1 += tableau(int).b[j] * F[j][k]
                y2 += tableau(int).b̂[j] * F[j][k]
            end
            b[ndims(int) * nstages(int) + k] = θ[k] - p̄[k] - timestep(int) * (y1 + y2)
        end
    end
end


# Compute stages of implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::IntegratorIRKimplicit) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::IntegratorIRKimplicit) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    update!(solstep(int), cache(int, DT).V, cache(int, DT).F, tableau(int), timestep(int))
end


function integrate_step!(int::IntegratorIRKimplicit)
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(nlsolution(int), int)

    # update one-form for next step
    # functions(problem).ϑ(solstep.p, solstep.t, solstep.q, solstep.v)
end

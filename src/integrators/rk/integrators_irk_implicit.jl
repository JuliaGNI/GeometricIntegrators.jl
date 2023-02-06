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

        q = zeros(DT,D)
        v = zeros(DT,D)
        θ = zeros(DT,D)
        f = zeros(DT,D)

        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        
        new(x, q, v, θ, f, Q, V, Θ, F)
    end
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


function initsolver(::Newton, solstep::SolutionStepPODE{DT}, problem::Union{IODEProblem,LODEProblem}, method::IRKMethod, caches::CacheDict) where {DT}
    # create wrapper function f!(b,x)
    f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches)

    # create nonlinear solver
    NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
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
        initialguess!(solstep.t̄[1] + timestep(problem) * tableau(method).c[i], cache.Q[i], cache.Θ[i], cache.V[i], cache.F[i], solstep, problem, iguess)
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(method)
        offset = ndims(problem)*(i-1)
        for k in 1:ndims(problem)
            cache.x[offset+k] = cache.Θ[i][k] - solstep.p̄[1][k]
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


function compute_stages!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::IRKMethod,
    caches::CacheDict) where {ST,DT,TT}

    # temporary variables
    local y1::ST
    local y2::ST
    local tᵢ::TT

    # get cache for internal stages
    cache = caches[ST]
    Q, V, Θ, F = cache.Q, cache.V, cache.Θ, cache.F

    # copy x to V
    for i in eachindex(V)
        for k in eachindex(V[i])
            V[i][k] = x[ndims(problem)*(i-1) + k]
        end
    end

    # copy x to q
    if implicit_update(method)
        for k in eachindex(cache.q)
            cache.q[k] = x[ndims(problem) * nstages(tableau(method)) + k]
        end
    end

    # compute Q = q + Δt A V
    for i in eachindex(Q)
        for k in eachindex(Q[i])
            y1 = y2 = 0
            for j in eachindex(V)
                y1 += tableau(method).a[i,j] * V[j][k]
                y2 += tableau(method).â[i,j] * V[j][k]
            end
            Q[i][k] = solstep.q̄[1][k] + timestep(problem) * (y1 + y2)
        end
    end

    # compute Θ = ϑ(Q) and F = f(Q,V)
    for i in eachindex(Θ,F)
        tᵢ = solstep.t̄[1] + timestep(problem) * tableau(method).c[i]
        functions(problem).ϑ(Θ[i], tᵢ, Q[i], V[i])
        functions(problem).f(F[i], tᵢ, Q[i], V[i])
    end

    # compute θ = ϑ(q)
    if implicit_update(method)
        functions(problem).ϑ(cache.θ, solstep.t, cache.q, cache.v)
    end
end

# Compute stages of implicit Runge-Kutta methods.
function function_stages!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::IRKMethod,
    caches::CacheDict) where {ST}

    # temporary variables
    local y1::ST
    local y2::ST

    # get cache
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b for internal stages
    for i in eachindex(cache.Θ)
        for k in eachindex(cache.Θ[i])
            y1 = y2 = 0
            for j in eachindex(cache.F)
                y1 += tableau(method).a[i,j] * cache.F[j][k]
                y2 += tableau(method).â[i,j] * cache.F[j][k]
            end
            b[ndims(problem)*(i-1) + k] = cache.Θ[i][k] - solstep.p̄[1][k] - timestep(problem) * (y1 + y2)
        end
    end

    # compute b for update
    if implicit_update(method)
        for k in eachindex(cache.θ)
            y1 = 0
            y2 = 0
            for j in eachindex(cache.F)
                y1 += tableau(method).b[j] * cache.F[j][k]
                y2 += tableau(method).b̂[j] * cache.F[j][k]
            end
            b[ndims(problem) * nstages(tableau(method)) + k] = cache.θ[k] - solstep.p̄[1][k] - timestep(problem) * (y1 + y2)
        end
    end
end


function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, 
    method::IRKMethod,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute vector field at internal stages
    compute_stages!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    update_solution!(solstep.q, solstep.q̄[1], solstep.q̃, caches[DT].V, tableau(method).b, tableau(method).b̂, timestep(problem))

    # update vector field for initial guess
    update_vector_fields!(solstep, problem)

    # update one-form for next step
    # functions(problem).ϑ(solstep.p, solstep.t, solstep.q, solstep.v)
end

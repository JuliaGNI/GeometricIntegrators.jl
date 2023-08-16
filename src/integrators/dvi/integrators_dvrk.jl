@doc raw"""
Degenerate variational Runge-Kutta integrator cache.

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
struct IntegratorCacheDVRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Θ::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function IntegratorCacheDVRK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        new(zeros(DT, D*(S+1)), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Θ, F)
    end
end

nlsolution(cache::IntegratorCacheDVRK) = cache.x

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::DVRK; kwargs...) where {ST}
    IntegratorCacheDVRK{ST, ndims(problem), nstages(tableau(method))}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::DVRK) = IntegratorCacheDVRK{ST, ndims(problem)}



@doc raw"""
Degenerate Variational Runge-Kutta integrator for noncanonical
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
"""
const IntegratorDVRK{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:DVRK}


default_solver(::DVRK) = Newton()
default_iguess(::DVRK) = HermiteExtrapolation()


function Base.show(io::IO, int::IntegratorDVRK)
    print(io, "\nRunge-Kutta Integrator for Degenerate Lagrangians with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(tableau(int)))
end


function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::Union{IODEProblem,LODEProblem},
    method::DVRK, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    # get cache and dimension
    cache = caches[DT]

    # compute initial guess for internal stages
    for i in eachstage(tableau(method))
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).c[i], cache.Q[i], cache.Θ[i], cache.V[i], cache.F[i], solstep, problem, iguess)
    end
    for i in eachstage(tableau(method))
        for k in 1:ndims(problem)
            cache.x[ndims(problem)*(i-1)+k] = cache.V[i][k]
        end
    end
    
    # compute initial guess for solution
    initialguess!(solstep.t, cache.q, cache.θ, cache.v, cache.f, solstep, problem, iguess)

    for k in 1:ndims(problem)
        cache.x[ndims(problem) * nstages(tableau(method)) + k] = cache.q[k]
    end
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::DVRK,
    caches::CacheDict) where {ST,DT,TT}

    # get cache and dimension
    cache = caches[ST]
    D = ndims(problem)
    S = nstages(tableau(method))

    # temporary variables
    local y1::ST
    local y2::ST
    local tᵢ::TT

    # copy x to V
    for i in eachindex(cache.V)
        for k in eachindex(cache.V[i])
            cache.V[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to q
    for k in eachindex(cache.q)
        cache.q[k] = x[D*S+k]
    end

    # compute Q = q + Δt A V, Θ = ϑ(Q), F = f(Q,V)
    for i in eachindex(cache.Q, cache.F, cache.Θ)
        tᵢ = solstep.t̄ + timestep(problem) * tableau(method).c[i]
        for k in eachindex(cache.Q[i])
            y1 = 0
            y2 = 0
            for j in eachindex(cache.V)
                y1 += tableau(method).a[i,j] * cache.V[j][k]
                y2 += tableau(method).â[i,j] * cache.V[j][k]
            end
            cache.Q[i][k] = solstep.q̄[k] + timestep(problem) * (y1 + y2)
        end
        functions(problem).ϑ(cache.Θ[i], tᵢ, cache.Q[i], cache.V[i])
        functions(problem).f(cache.F[i], tᵢ, cache.Q[i], cache.V[i])
    end

    # compute q̄ = q + Δt B V, Θ = ϑ(q̄)
    functions(problem).ϑ(cache.θ, solstep.t, cache.q, cache.v)
end


# Compute stages of fully implicit Runge-Kutta methods.
function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::DVRK,
    caches::CacheDict) where {ST}

    # get cache
    cache = caches[ST]
    D = ndims(problem)
    S = nstages(tableau(method))

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)
    
    # compute b
    for i in eachindex(cache.Θ)
        for k in eachindex(cache.Θ[i])
            y1 = 0
            y2 = 0
            for j in eachindex(cache.F)
                y1 += tableau(method).a[i,j] * cache.F[j][k]
                y2 += tableau(method).â[i,j] * cache.F[j][k]
            end
            b[D*(i-1)+k] = cache.Θ[i][k] - solstep.p̄[k] - timestep(problem) * (y1 + y2)
        end
    end
    for k in 1:div(D,2)
        y1 = 0
        y2 = 0
        for j in eachindex(cache.F)
            y1 += tableau(method).b[j] * cache.F[j][k]
            y2 += tableau(method).b̂[j] * cache.F[j][k]
        end
        b[D*S+k] = cache.θ[k] - solstep.p̄[k] - timestep(problem) * (y1 + y2)
    end
    for k in 1:div(D,2)
        y1 = 0
        y2 = 0
        for j in eachindex(cache.V)
            y1 += tableau(method).b[j] * cache.V[j][k]
            y2 += tableau(method).b̂[j] * cache.V[j][k]
        end
        b[D*S+div(D,2)+k] = cache.q[k] - solstep.q̄[k] - timestep(problem) * (y1 + y2)
    end
end


function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem{DT,TT},LODEProblem{DT,TT}},
    method::DVRK,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, (b,x) -> residual!(b, x, solstep, problem, method, caches), solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    components!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    solstep.q .= caches[DT].q
    solstep.p .= caches[DT].θ
end

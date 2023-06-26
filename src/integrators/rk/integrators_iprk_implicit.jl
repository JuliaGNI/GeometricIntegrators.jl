@doc raw"""
Implicit partitioned Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of q
* `p̃`: initial guess of p
* `ṽ`: initial guess of v
* `f̃`: initial guess of f
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of q
* `P`: internal stages of p
* `V`: internal stages of v
* `F`: internal stages of f
* `Y`: vector field of internal stages of q
* `Z`: vector field of internal stages of p
"""
struct IntegratorCacheIPRKimplicit{ST,D,S} <: PODEIntegratorCache{ST,D}
    x::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    function IntegratorCacheIPRKimplicit{ST,D,S}() where {ST,D,S}
        # create solver vector
        x = zeros(ST, D*S)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)

        new(x, Q, P, V, F)
    end
end

nlsolution(cache::IntegratorCacheIPRKimplicit) = cache.x


function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::IPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheIPRKimplicit{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::IPRK) = IntegratorCacheIPRKimplicit{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Implicit partitioned Runge-Kutta integrator solving the system
```math
\begin{aligned}
P_{n,i} &= \vartheta (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
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
const IntegratorIPRKimplicit{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:IPRK}

solversize(problem::Union{IODEProblem,LODEProblem}, method::IPRK) =
    ndims(problem) * nstages(method)


function initsolver(::Newton, solstep::SolutionStepPODE{DT}, problem::Union{IODEProblem,LODEProblem}, method::IPRKMethod, caches::CacheDict) where {DT}
    # create wrapper function f!(b,x)
    f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches)

    # create nonlinear solver
    NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end


function Base.show(io::IO, int::IntegratorIPRKimplicit)
    print(io, "\nPartitioned Runge-Kutta Integrator for Implicit Equations with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end


function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::Union{IODEProblem,LODEProblem},
    method::IPRKMethod, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    local x = caches[DT].x
    local Q = caches[DT].Q
    local P = caches[DT].P
    local V = caches[DT].V
    local F = caches[DT].F

    # compute initial guess for internal stages
    for i in eachstage(method)
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).p.c[i], Q[i], P[i], V[i], F[i], solstep, problem, iguess)
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(method)
        offset = ndims(problem)*(i-1)
        for k in 1:ndims(problem)
            x[offset+k] = P[i][k] - solstep.p̄[k]
            for j in eachstage(method)
                x[offset+k] -= timestep(problem) * tableau(method).p.a[i,j] * F[j][k]
            end
        end
    end
end


function compute_stages!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::IPRKMethod,
    caches::CacheDict) where {ST,DT,TT}

    # get cache for internal stages
    local Q = caches[ST].Q
    local P = caches[ST].P
    local V = caches[ST].V
    local F = caches[ST].F
    local D = ndims(problem)

    # temporary variables
    local tqᵢ::TT
    local tpᵢ::TT
    local y1::ST
    local y2::ST

    for i in eachstage(method)
        for k in 1:D
            # copy y to V
            V[i][k] = x[D*(i-1)+k]

            # compute Q
            y1 = y2 = 0
            for j in eachstage(method)
                y1 += tableau(method).q.a[i,j] * V[j][k]
                y2 += tableau(method).q.â[i,j] * V[j][k]
            end
            Q[i][k] = solstep.q̄[k] + timestep(problem) * (y1 + y2)
        end

        # compute time of internal stage
        tqᵢ = solstep.t̄ + timestep(problem) * tableau(method).q.c[i]
        tpᵢ = solstep.t̄ + timestep(problem) * tableau(method).p.c[i]

        # compute ϑ(Q,V) and f(Q,V)
        functions(problem).ϑ(P[i], tqᵢ, Q[i], V[i])
        functions(problem).f(F[i], tpᵢ, Q[i], V[i])
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function function_stages!(
    b::AbstractVector{ST},
    x::AbstractVector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::IPRKMethod,
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local P = caches[ST].P
    local F = caches[ST].F
    local D = ndims(problem)

    # temporary variables
    local z1::ST
    local z2::ST

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b = - [(Y-AV), (Z-AF)]
    for i in eachstage(method)
        for k in 1:D
            z1 = z2 = 0
            for j in eachstage(method)
                z1 += tableau(method).p.a[i,j] * F[j][k]
                z2 += tableau(method).p.â[i,j] * F[j][k]
            end
            b[D*(i-1)+k] = P[i][k] - solstep.p̄[k] - timestep(problem) * (z1 + z2)
        end
    end
end


function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem{DT,TT},LODEProblem{DT,TT}},
    method::IPRK,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, (b,x) -> function_stages!(b, x, solstep, problem, method, caches), solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute vector fields at internal stages
    compute_stages!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    update!(solstep, caches[DT].V, caches[DT].F, tableau(method), timestep(problem))
end

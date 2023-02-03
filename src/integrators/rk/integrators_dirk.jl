"""
Diagonally implicit Runge-Kutta integrator cache.
"""
struct IntegratorCacheDIRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    x::Vector{Vector{DT}}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCacheDIRK{DT,D,S}() where {DT,D,S}
        x = create_internal_stage_vector(DT, D, S)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        new(x, Q, V, Y)
    end
end

function Cache{ST}(problem::GeometricProblem, method::DIRKMethod; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheDIRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::GeometricProblem, method::DIRKMethod) = IntegratorCacheDIRK{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Diagonally implicit Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
"""
const IntegratorDIRK{DT,TT} = Integrator{<:ODEProblem{DT,TT}, <:DIRKMethod}

initmethod(method::DIRK) = method
initmethod(method::DIRKMethod) = DIRK(method)

default_solver(::DIRKMethod) = Newton()
default_iguess(::DIRKMethod) = HermiteExtrapolation()


function Base.show(io::IO, int::IntegratorDIRK)
    print(io, "\nDiagonally Implicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


struct SingleStageSolvers{ST} <: NonlinearSolver
    solvers::ST
    SingleStageSolvers(solvers...) = new{typeof(solvers)}(solvers)
end

Base.getindex(s::SingleStageSolvers, args...) = getindex(s.solvers, args...)


function initsolver(::Newton, solstep::SolutionStepODE{DT}, problem::ODEProblem, method::DIRKMethod, caches::CacheDict, i) where {DT}
    # create wrapper function f!(b,x)
    f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches, i)

    # create nonlinear solver
    NewtonSolver(zero(caches[DT].x[i]), zero(caches[DT].x[i]), f!; linesearch = Backtracking(), config = Options(min_iterations = 1))
end

function initsolver(solvermethod::Newton, solstep::SolutionStepODE, problem::ODEProblem, method::DIRKMethod, caches::CacheDict)
    SingleStageSolvers([initsolver(solvermethod, solstep, problem, method, caches, i) for i in 1:nstages(tableau(method))]...)
end


function initial_guess!(
    solstep::SolutionStepODE{DT}, 
    problem::ODEProblem,
    method::DIRKMethod, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}
    
    # obtain cache
    cache = caches[DT]

    for i in eachstage(method)
        initialguess!(solstep.t̄[1] + timestep(problem) * tableau(method).c[i], cache.Q[i], cache.V[i], solstep, problem, iguess)
    end

    for i in eachstage(method)
        for k in 1:ndims(problem)
            cache.x[i][k] = 0
            for j in eachstage(method)
                cache.x[i][k] += timestep(problem) * tableau(method).a[i,j] * cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, solstep::SolutionStepODE{DT,TT}, problem::ODEProblem, method::DIRKMethod, caches::CacheDict, i) where {ST,DT,TT}
    local tᵢ::TT

    # get cache for internal stages
    local Q::Vector{Vector{ST}} = caches[ST].Q
    local V::Vector{Vector{ST}} = caches[ST].V
    local Y::Vector{Vector{ST}} = caches[ST].Y

    # copy x to Y and compute Q = q + Y
    for k in 1:ndims(problem)
        Y[i][k] = x[k]
        Q[i][k] = solstep.q̄[1][k] + Y[i][k]
    end

    # compute V = v(Q)
    tᵢ = solstep.t̄[1] + timestep(problem) * tableau(method).c[i]
    functions(problem).v(V[i], tᵢ, Q[i])
end


# Compute stages of diagonally implicit Runge-Kutta methods.
function function_stages!(b::Vector{ST}, x::Vector{ST}, solstep::SolutionStepODE, problem::ODEProblem, method::DIRKMethod, caches::CacheDict, i) where {ST}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    local V::Vector{Vector{ST}} = caches[ST].V
    local Y::Vector{Vector{ST}} = caches[ST].Y

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches, i)

    # compute b = - (Y-AV)
    for k in 1:ndims(problem)
        y1 = tableau(method).a[i,i] * V[i][k]
        y2 = tableau(method).â[i,i] * V[i][k]
        for j in 1:i-1
            y1 += tableau(method).a[i,j] * V[j][k]
            y2 += tableau(method).â[i,j] * V[j][k]
        end
        b[k] = - Y[i][k] + timestep(problem) * (y1 + y2)
    end
end


function integrate_step!(
    solstep::SolutionStepODE{DT,TT},
    problem::ODEProblem{DT,TT}, 
    method::DIRKMethod,
    caches::CacheDict,
    solvers::SingleStageSolvers) where {DT,TT}

    # consecutively solve for all stages
    for i in eachstage(method)
        # call nonlinear solver
        solve!(caches[DT].x[i], solvers[i])

        # print solver status
        # println(status(solvers[i]))

        # check if solution contains NaNs or error bounds are violated
        # println(meets_stopping_criteria(status(solvers[i])))

        # compute vector field at internal stages
        compute_stages!(caches[DT].x[i], solstep, problem, method, caches, i)
    end

    # compute final update
    update_solution!(solstep.q, solstep.q̄[1], solstep.q̃, caches[DT].V, tableau(method).b, tableau(method).b̂, timestep(problem))

    # update vector field for initial guess
    update_vector_fields!(solstep, problem)
end

@doc raw"""
Degenerate variational integrator cache.

### Fields

* `q`: internal stages of solution
* `v`: internal stages of vector field
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
"""
struct IntegratorCacheDVIA{DT,D} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    v̄::Vector{DT}
    θ̄::Vector{DT}
    f̄::Vector{DT}

    function IntegratorCacheDVIA{DT,D}() where {DT,D}
        new(zeros(DT,2D), 
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
                         zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::DVIA; kwargs...) where {ST}
    IntegratorCacheDVIA{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::DVIA) = IntegratorCacheDVIA{ST, ndims(problem)}


"""
Symplectic Euler-A Degenerate Variational Integrator.
"""
const IntegratorDVIA{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:DVIA}


default_solver(::DVIA) = Newton()
default_iguess(::DVIA) = HermiteExtrapolation()


function initsolver(::Newton, solstep::SolutionStepPODE{DT}, problem::Union{IODEProblem,LODEProblem}, method::DVIA, caches::CacheDict) where {DT}
    # create wrapper function f!(b,x)
    f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches)

    # create nonlinear solver
    NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end


function Base.show(io::IO, int::IntegratorDVIA)
    print(io, "\nDegenerate Variational Integrator (Euler-A) with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::Union{IODEProblem,LODEProblem},
    method::DVIA, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    # get cache and dimension
    cache = caches[DT]
    D = ndims(problem)

    # compute initial guess for solution
    initialguess!(solstep.t, cache.q, cache.θ, cache.v, cache.f, solstep, problem, iguess)

    cache.x[1:D] .= cache.q

    for k in 1:div(D,2)
        cache.x[D+k] = cache.v[k]
        cache.x[D+div(D,2)+k] = solstep.v[k]
    end
end


function compute_stages!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::DVIA,
    caches::CacheDict) where {ST,DT,TT}

    # get cache and dimension
    cache = caches[ST]
    D = ndims(problem)

    # copy x to q
    cache.q .= x[1:D]

    # copy x to v and v̄
    for k in 1:div(D,2)
        cache.v[k] = x[D+k]
        cache.v̄[k] = x[D+div(D,2)+k]
        cache.v[div(D,2)+k] = 0
        cache.v̄[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    functions(problem).f(cache.f, solstep.t, cache.q, cache.v)
    functions(problem).f(cache.f̄, solstep.t̄[1], solstep.q̄[1], cache.v̄)
 
    # compute Θ = ϑ(q,v)
    functions(problem).ϑ(cache.θ, solstep.t, cache.q, cache.v)
    # functions(problem).ϑ(cache.θ̄, solstep.t̄[1], solstep.q̄[1], cache.v̄)
end


function function_stages!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::DVIA,
    caches::CacheDict) where {ST}

    # get cache
    cache = caches[ST]
    D = ndims(problem)

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b
    for k in 1:div(D,2)
        b[k]   = cache.θ[k] - solstep.p̄[1][k] - timestep(problem) * cache.f̄[k]
        b[D+k] = cache.q[k] - solstep.q̄[1][k] - timestep(problem) * cache.v[k]
    end

    for k in div(D,2)+1:D
        b[k]   = timestep(problem) * cache.f̄[k]
        b[D+k] = timestep(problem) * cache.f[k]
    end
end


function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem{DT,TT},LODEProblem{DT,TT}},
    method::DVIA,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, solver)

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    solstep.q .= caches[DT].q
    solstep.v .= caches[DT].v
    solstep.p .= caches[DT].θ
    solstep.f .= caches[DT].f

    # update vector field for initial guess
    update_vector_fields!(solstep, problem)
end

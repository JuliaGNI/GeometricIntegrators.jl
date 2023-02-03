@doc raw"""
Degenerate variational integrator cache.

### Fields

* `q`: internal stages of solution
* `v`: internal stages of vector field
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
"""
struct IntegratorCacheCTDVI{DT,D} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    q⁻::Vector{DT}
    θ⁻::Vector{DT}
    f⁻::Vector{DT}

    q⁺::Vector{DT}
    θ⁺::Vector{DT}
    f⁺::Vector{DT}

    function IntegratorCacheCTDVI{DT,D}() where {DT,D}
        new(zeros(DT,2D), 
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::CTDVI; kwargs...) where {ST}
    IntegratorCacheCTDVI{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::CTDVI) = IntegratorCacheCTDVI{ST, ndims(problem)}


"""
Trapezoidal Degenerate Variational Integrator.
"""
const IntegratorCTDVI{DT,TT} = Integrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:CTDVI}


default_solver(::CTDVI) = Newton()
default_iguess(::CTDVI) = HermiteExtrapolation()


function initsolver(::Newton, solstep::SolutionStepPODE{DT}, problem::Union{IODEProblem,LODEProblem}, method::CTDVI, caches::CacheDict) where {DT}
    # create wrapper function f!(b,x)
    f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches)

    # create nonlinear solver
    NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1, x_abstol = 8eps(), f_abstol = 8eps()))
end


function Base.show(io::IO, int::IntegratorCTDVI)
    print(io, "\nTrapezoidal Degenerate Variational Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end



function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::Union{IODEProblem,LODEProblem},
    method::CTDVI, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    # get cache and dimension
    cache = caches[DT]
    D = ndims(problem)

    # compute initial guess for solution q(n+1)
    initialguess!(solstep.t, cache.q, cache.θ, cache.v, cache.f, solstep, problem, iguess)

    cache.x[1:D] .= cache.q

    # compute initial guess for solution q(n+1/2)
    initialguess!((solstep.t + solstep.t̄[1])/2, cache.q, cache.θ, cache.v, cache.f, solstep, problem, iguess)

    offset_v = D
    offset_x = D + div(D,2)
    for k in 1:div(D,2)
        cache.x[offset_v+k] = cache.v[k]              # v¹(n+1/2)
        cache.x[offset_x+k] = cache.q[div(D,2)+k]     # q²(n+1/2)
    end
end


function compute_stages!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::CTDVI,
    caches::CacheDict) where {ST,DT,TT}

    # get cache and dimension
    cache = caches[ST]
    D = ndims(problem)

    # set some local variables for convenience and clarity
    local t⁻ = solstep.t̄[1]
    local t⁺ = solstep.t̄[1] + timestep(problem)
    
    # copy x to q
    cache.q .= x[1:D]

    # copy x to q⁻, q⁺ and v
    for k in 1:div(D,2)
        cache.q⁻[k] = solstep.q̄[1][k]
        cache.q⁺[k] = cache.q[k]

        cache.q⁻[div(D,2)+k] = x[D+div(D,2)+k]
        cache.q⁺[div(D,2)+k] = x[D+div(D,2)+k]

        cache.v[k]          = x[D+k]
        cache.v[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    functions(problem).f(cache.f⁻, t⁻, cache.q⁻, cache.v)
    functions(problem).f(cache.f⁺, t⁺, cache.q⁺, cache.v)
 
    # compute Θ = ϑ(q,v)
    functions(problem).ϑ(cache.θ⁻, t⁻, cache.q⁻, cache.v)
    functions(problem).ϑ(cache.θ⁺, t⁺, cache.q⁺, cache.v)
    functions(problem).ϑ(cache.θ,  t⁺, cache.q, cache.v)
end


function function_stages!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::CTDVI,
    caches::CacheDict) where {ST}

    # get cache for internal stages
    cache = caches[ST]
    D = ndims(problem)

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b
    b[1:D] .= (cache.θ⁻ .+ cache.θ⁺) ./ 2 .- solstep.p̄[1] .- timestep(problem) .* cache.f⁻ ./ 2
    
    for k in 1:div(D,2)
        b[D+k]          = cache.q[k] - solstep.q̄[1][k] - timestep(problem) * cache.v[k]
        b[D+div(D,2)+k] = cache.θ[k] - solstep.p̄[1][k] - timestep(problem) * (cache.f⁻[k] + cache.f⁺[k]) / 2
    end
end


function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem{DT,TT},LODEProblem{DT,TT}},
    method::CTDVI,
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
    solstep.p .= caches[DT].θ

    # update vector field for initial guess
    update_vector_fields!(solstep, problem)
end

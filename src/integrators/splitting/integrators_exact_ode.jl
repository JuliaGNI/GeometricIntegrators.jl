
"Exact solution of an ODE."
struct ExactSolutionODE{T, QT <: Base.Callable} <: ODEMethod
    c::T
    q::QT
end

function Cache{ST}(problem::SODEProblem, method::ExactSolutionODE; kwargs...) where {ST}
    IntegratorCacheSplitting{ST, typeof(timestep(problem)), ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::SODEProblem, ::ExactSolutionODE) = IntegratorCacheSplitting{ST, typeof(timestep(problem)), ndims(problem)}


function integrate_step!(
    solstep::SolutionStepODE{DT,TT},
    problem::SODEProblem{DT,TT},
    method::ExactSolutionODE,
    caches::CacheDict,
    ::NoSolver) where {DT,TT}
    
    # copy previous solution
    caches[DT].q .= solstep.q

    # compute new solution
    method.q(solstep.q, sol.t̄[1] + timestep(problem) * method.c, caches[DT].q, solstep.t̄[1])
end

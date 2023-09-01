"""
Implements the [`GeometricIntegrator`](@ref) interface for the exact solution of a substep in [`SODE`](@ref).
"""
struct ExactSolution <: GeometricMethod end

function Cache{ST}(problem::SubstepProblem, method::ExactSolution; kwargs...) where {ST}
    SplittingCache{ST, typeof(timestep(problem)), ndims(problem)}(initial_conditions(problem).q; kwargs...)
end

@inline CacheType(ST, problem::SubstepProblem, ::ExactSolution) = SplittingCache{ST, typeof(timestep(problem)), ndims(problem)}


function integrate_step!(
    solstep::SolutionStepODE{DT,TT},
    problem::SubstepProblem,
    method::ExactSolution,
    caches::CacheDict,
    ::NoSolver) where {DT,TT}
    
    # copy previous solution
    caches[DT].q .= solstep.q

    # compute new solution
    solutions(problem).q(solstep.q, solstep.t̄ + timestep(problem), caches[DT].q, solstep.t̄)
end

"""
Implements the [`GeometricIntegrator`](@ref) interface for the exact solution of a substep in [`SODE`](@ref).
"""
struct ExactSolution <: GeometricMethod end

function Cache{ST}(problem::SubstepProblem, method::ExactSolution; kwargs...) where {ST}
    SplittingCache{ST, typeof(timestep(problem)), ndims(problem)}(initial_conditions(problem).q; kwargs...)
end

@inline CacheType(ST, problem::SubstepProblem, ::ExactSolution) = SplittingCache{ST, typeof(timestep(problem)), ndims(problem)}


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ExactSolution, <:SubstepProblem})
    # copy previous solution
    cache(int).q .= sol.q

    # compute new solution
    solutions(problem(int)).q(sol.q, history.t[1] + timestep(int), cache(int).q, history.t[1], params)
end

function integrate_step!(int::GeometricIntegrator{<:ExactSolution, <:SubstepProblem})
    integrate_step!(current(solstep(int)), history(solstep(int)), parameters(solstep(int)), int)
end

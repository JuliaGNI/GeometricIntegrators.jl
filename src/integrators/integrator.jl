#*****************************************************************************#
# Legacy code (we need to see which parts of this are still needed)           #
#*****************************************************************************#

# Cache{DT}(int::Integrator) where {DT} = Cache{DT}(int.problem, int.method)

# cache(::GeometricIntegrator, DT) = nothing
# cache(::GeometricIntegrator) = nothing

# initialize!(::SolutionStep, ::AbstractProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
# initialize!(int::DeterministicIntegrator) = initialize!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))

# function residual!(b::AbstractVector, x::AbstractVector, int::DeterministicIntegrator)
#     residual!(b, x, current(solstep(int)), history(solstep(int)), parameters(solstep(int)), int)
# end

# components!(x::AbstractVector, int::DeterministicIntegrator) = components!(x, solstep(int), problem(int), method(int), caches(int))

# GeometricBase.parameters(integrator::GeometricIntegrator) = error("parameters() not implemented for ", typeof(integrator))
# GeometricBase.equations(integrator::GeometricIntegrator) = error("equations() not implemented for ", typeof(integrator))
# GeometricBase.equation(integrator::GeometricIntegrator) = error("equation() not implemented for ", typeof(integrator))
# GeometricBase.equation(integrator::GeometricIntegrator, i::Int) = error("equation() not implemented for ", typeof(integrator))
# GeometricBase.timestep(integrator::GeometricIntegrator) = error("timestep() not implemented for ", typeof(integrator))
# Base.ndims(integrator::GeometricIntegrator) = error("ndims() not implemented for ", typeof(integrator))
# GeometricBase.nconstraints(integrator::GeometricIntegrator) = error("nconstraints() not implemented for ", typeof(integrator))
# # nstages(integrator::GeometricIntegrator) = error("nstages() not implemented for ", typeof(integrator))

# eachdim(integrator::GeometricIntegrator) = 1:ndims(integrator)


# Create SolutionStep with internal variables of integrator.
# function SolutionStep(solution::GeometricSolution, integrator::GeometricIntegrator)
#     SolutionStep(solution, internal_variables(integrator))
# end


# abstract type Parameters{DT,TT} end

# residual!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("residual!() not implemented for ", PT)
# solution_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("solution_stages!() not implemented for ", PT)

# initialize!(::GeometricIntegrator, ::SolutionStep) = nothing

# cache(::GeometricIntegrator) = missing

# """
# Performs one time step with a given integrator.

# ```julia
# integrate_step!(integrator::Integrator, solstep::SolutionStep)
# ```

# The function accepts two arguments: an integrator and an appropriate [`SolutionStep`](@ref),
# which contains the state of the system at the beginning and the end of the time step and possibly
# additional information like solver output or the solution at internal stages of a Runge-Kutta
# method.
# """
# integrate_step!(integrator::GeometricIntegrator, ::SolutionStep) = error("integrate_step!() not implemented for ", typeof(integrator))

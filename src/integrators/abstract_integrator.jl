
abstract type DeterministicIntegrator <: AbstractIntegrator end


struct NoSolver <: SolverMethod end

default_options() = Options(
    min_iterations = 1,
    x_abstol = 8eps(),
    f_abstol = 8eps(),
)

default_solver(::GeometricMethod) = NoSolver()
# default_solver(::Nothing) = NoSolver()

default_iguess(::GeometricMethod) = NoInitialGuess()
# default_iguess(::Nothing) = NoInitialGuess()

# default_projection(::GeometricMethod) = NoProjection()

initmethod(method::GeometricMethod) = method
initmethod(method::GeometricMethod, ::AbstractProblem) = initmethod(method)
initsolver(::SolverMethod, ::GeometricMethod, ::CacheDict) = NoSolver()

# create nonlinear solver
function initsolver(::NewtonMethod, ::GeometricMethod, caches::CacheDict; config = default_options(), kwargs...)
    x = zero(nlsolution(caches))
    y = zero(nlsolution(caches))
    NewtonSolver(x, y; linesearch = Backtracking(), config = config, kwargs...)
end




#*****************************************************************************#
# General integration functions for all integrators                           #
#*****************************************************************************#


# Apply integrator for ntime time steps and return solution.
function integrate(problem::AbstractProblem, method::GeometricMethod; kwargs...)
    integrator = GeometricIntegrator(problem, method; kwargs...)
    solution = Solution(problem; kwargs...)
    integrate!(solution, integrator)
end

function integrate(problems::GeometricEnsemble, method::GeometricMethod; kwargs...)
    solutions = Solution(problems; kwargs...)

    for (problem, solution) in zip(problems, solutions)
        integrator = GeometricIntegrator(problem, method; kwargs...)
        integrate!(solution, integrator)
    end

    return solutions
end

integrate_step!(int::DeterministicIntegrator) = integrate_step!(solstep(int), problem(int), method(int), caches(int), solver(int))


# Parts of one integration step that are common to deterministic and stochastic equations.
# function integrate!(solstep::SolutionStep, problem::EquationProblem, method::GeometricMethod, caches::CacheDict, solver::Union{SolverMethod, NonlinearSolver}, iguess::Union{InitialGuess,Extrapolation})
function integrate!(int::DeterministicIntegrator)
    # reset solution step
    reset!(solstep(int), timestep(int))

    # compute initial guess
    initial_guess!(int)

    # integrate one initial condition for one time step
    integrate_step!(int)

    # take care of periodic solutions
    cut_periodic_solution!(solstep(int), periodicity(problem(int)))

    # update vector field for initial guess
    Solutions.update_vector_fields!(solstep(int), problem(int))

    return solstep(int)
end

# Integrate equation for time steps n with n₁ ≤ n ≤ n₂.
"""
Solve for time steps n with n₁ ≤ n ≤ n₂.
```julia
integrate!(solution, integrator, n₁, n₂)
```
"""
function integrate!(sol::GeometricSolution, int::DeterministicIntegrator, n₁::Int, n₂::Int)
    # check time steps range for consistency
    @assert n₁ ≥ 1
    @assert n₂ ≥ n₁
    @assert n₂ ≤ ntime(sol)

    # copy initial condition from solution and initialize
    copy!(solstep(int), sol[n₁-1])
    initialize!(int)
    Solutions.initialize!(solstep(int), problem(int), Solutions.default_extrapolation())

    # loop over time steps
    for n in n₁:n₂
        # integrate one step and copy solution from cache to solution
        sol[n] = integrate!(int)

        # try
        #     sol[n] = integrate!(int)
        # catch ex
        #     tstr = " in time step " * string(n)
        #
        #     if m₁ ≠ m₂
        #         tstr *= " for initial condition " * string(m)
        #     end
        #
        #     tstr *= "."
        #
        #     if isa(ex, DomainError)
        #         @warn("Domain error" * tstr)
        #     elseif isa(ex, ErrorException)
        #         @warn("Simulation exited early" * tstr)
        #         @warn(ex.msg)
        #     else
        #         @warn(string(typeof(ex)) * tstr)
        #         throw(ex)
        #     end
        # end
    end

    return sol
end

"""
Solve for all time steps n:
```julia
integrate!(solution, integrator)
```
"""
function integrate!(sol::GeometricSolution, int::DeterministicIntegrator)
    integrate!(sol, int, 1, ntime(sol))
end



#*****************************************************************************#
# Legacy code (we need to see which parts of this are still needed)           #
#*****************************************************************************#



# abstract type GeometricIntegrator{dType, tType} end

# abstract type DeterministicIntegrator{dType, tType} <: GeometricIntegrator{dType, tType} end

# GeometricBase.parameters(integrator::GeometricIntegrator) = error("parameters() not implemented for ", typeof(integrator))
# GeometricBase.equations(integrator::GeometricIntegrator) = error("equations() not implemented for ", typeof(integrator))
# GeometricBase.equation(integrator::GeometricIntegrator) = error("equation() not implemented for ", typeof(integrator))
# GeometricBase.equation(integrator::GeometricIntegrator, i::Int) = error("equation() not implemented for ", typeof(integrator))
# GeometricBase.timestep(integrator::GeometricIntegrator) = error("timestep() not implemented for ", typeof(integrator))
# Base.ndims(integrator::GeometricIntegrator) = error("ndims() not implemented for ", typeof(integrator))
# GeometricBase.nconstraints(integrator::GeometricIntegrator) = error("nconstraints() not implemented for ", typeof(integrator))
# # nstages(integrator::GeometricIntegrator) = error("nstages() not implemented for ", typeof(integrator))

# eachdim(integrator::GeometricIntegrator) = 1:ndims(integrator)

"""
```julia
get_internal_variables(::Integrator) = NamedTuple()
```
Returns a `NamedTuple` containing all internal variables of an integrator that
shall be stored in an [`SolutionStep`](@ref). If there is no method for a
specific integrator implemented an empty `NamedTuple()` is returned.
"""
get_internal_variables(::DeterministicIntegrator) = NamedTuple()
get_internal_variables(::Nothing) = NamedTuple()


# Create SolutionStep with internal variables of integrator.
# function Solutions.SolutionStep(solution::GeometricSolution, integrator::GeometricIntegrator)
#     SolutionStep(solution, get_internal_variables(integrator))
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

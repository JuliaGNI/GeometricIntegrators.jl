
abstract type AbstractIntegrator end

abstract type DeterministicIntegrator <: AbstractIntegrator end

abstract type ODEIntegrator{dType, tType} <: DeterministicIntegrator end
abstract type DAEIntegrator{dType, tType} <: DeterministicIntegrator end
abstract type PODEIntegrator{dType, tType} <: DeterministicIntegrator end
abstract type PDAEIntegrator{dType, tType} <: DeterministicIntegrator end

abstract type IODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type IDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type HODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type HDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type LODEIntegrator{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type LDAEIntegrator{dType, tType} <: IDAEIntegrator{dType, tType} end



struct NoSolver <: SolverMethod end

default_solver(::GeometricMethod) = NoSolver()
# default_solver(::Nothing) = NoSolver()

default_iguess(::GeometricMethod) = NoInitialGuess()
# default_iguess(::Nothing) = NoInitialGuess()

# default_projection(::GeometricMethod) = NoProjection()

initmethod(method::GeometricMethod) = method
initmethod(method::GeometricMethod, ::AbstractProblem) = initmethod(method)
initsolver(::SolverMethod, ::GeometricMethod, ::CacheDict) = NoSolver()

# create nonlinear solver
function initsolver(::NewtonMethod, ::GeometricMethod, caches::CacheDict; kwargs...)
    x = zero(nlsolution(caches))
    y = zero(nlsolution(caches))
    NewtonSolver(x, y; linesearch = Backtracking(), config = Options(min_iterations = 1), kwargs...)
end




#*****************************************************************************#
# General integration functions for all integrators                           #
#*****************************************************************************#

"""
```julia
integrate(equation, integrator, ntime; kwargs...)
integrate(equation, tableau, Δt, ntime; kwargs...)
```

Integrate an `equation` with `integrator` for `ntime` time steps
and return the solution.
If a `tableau` and a time step `Δt` is passed instead of an
`integrator`, the appropriate integrator is created automatically.

Some convenience methods exist for the integration of ODEs,
```julia
integrate(v, x₀, tableau, Δt, ntime; t₀=0., kwargs...)
```
where `v` is the function for the vector field, `x₀` the initial condition
and `t₀` the initial time, and for PODEs
```julia
integrate(v, f, q₀, p₀, tableau, Δt, ntime; t₀=0., kwargs...)
```
with vector fields `v` and `f` and initial conditions `q₀` and `p₀`.
"""
function integrate end

# Apply integrator for ntime time steps and return solution.
function integrate(problem::AbstractProblem, method::GeometricMethod; kwargs...)
    integrator = Integrator(problem, method; kwargs...)
    solution = Solution(problem; kwargs...)
    integrate!(solution, integrator)
end

integrate_step!(int::AbstractIntegrator) = integrate_step!(solstep(int), problem(int), method(int), caches(int), solver(int))

"""
Solve for all time steps n:
```julia
integrate!(solution, integrator)
```
    
Solve for time steps n with n₁ ≤ n ≤ n₂.
```julia
integrate!(solution, integrator, n₁, n₂)
```

Solve one time step:
```julia
integrate!(integrator)
integrate!(::SolutionStep, ::GeometricProblem, ::GeometricMethod, ::CacheDict, ::SolverMethod, ::InitialGuess)
```
"""
function integrate! end

# Parts of one integration step that are common to deterministic and stochastic equations.
# function integrate!(solstep::SolutionStep, problem::GeometricProblem, method::GeometricMethod, caches::CacheDict, solver::Union{SolverMethod, NonlinearSolver}, iguess::Union{InitialGuess,Extrapolation})
function integrate!(int::AbstractIntegrator)
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
function integrate!(sol::GeometricSolution, int::AbstractIntegrator, n₁::Int, n₂::Int)
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

# Integrate `equation` for all time steps.
function integrate!(sol::GeometricSolution, int::AbstractIntegrator)
    integrate!(sol, int, 1, ntime(sol))
end




#*****************************************************************************#
# Legacy code (we need to see which parts of this are still needed)           #
#*****************************************************************************#



# abstract type AbstractIntegrator{dType, tType} end

# abstract type DeterministicIntegrator{dType, tType} <: AbstractIntegrator{dType, tType} end

# GeometricBase.parameters(integrator::AbstractIntegrator) = error("parameters() not implemented for ", typeof(integrator))
# GeometricBase.equations(integrator::AbstractIntegrator) = error("equations() not implemented for ", typeof(integrator))
# GeometricBase.equation(integrator::AbstractIntegrator) = error("equation() not implemented for ", typeof(integrator))
# GeometricBase.equation(integrator::AbstractIntegrator, i::Int) = error("equation() not implemented for ", typeof(integrator))
# GeometricBase.timestep(integrator::AbstractIntegrator) = error("timestep() not implemented for ", typeof(integrator))
# Base.ndims(integrator::AbstractIntegrator) = error("ndims() not implemented for ", typeof(integrator))
# GeometricBase.nconstraints(integrator::AbstractIntegrator) = error("nconstraints() not implemented for ", typeof(integrator))
# # nstages(integrator::AbstractIntegrator) = error("nstages() not implemented for ", typeof(integrator))

# eachdim(integrator::AbstractIntegrator) = 1:ndims(integrator)

"""
```julia
get_internal_variables(::Integrator) = NamedTuple()
```
Returns a `NamedTuple` containing all internal variables of an integrator that
shall be stored in an [`SolutionStep`](@ref). If there is no method for a
specific integrator implemented an empty `NamedTuple()` is returned.
"""
get_internal_variables(::AbstractIntegrator) = NamedTuple()
get_internal_variables(::Nothing) = NamedTuple()


# Create SolutionStep with internal variables of integrator.
# function Solutions.SolutionStep(solution::GeometricSolution, integrator::AbstractIntegrator)
#     SolutionStep(solution, get_internal_variables(integrator))
# end


# abstract type Parameters{DT,TT} end

# residual!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("residual!() not implemented for ", PT)
# solution_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("solution_stages!() not implemented for ", PT)

# initialize!(::AbstractIntegrator, ::SolutionStep) = nothing

# cache(::AbstractIntegrator) = missing

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
# integrate_step!(integrator::AbstractIntegrator, ::SolutionStep) = error("integrate_step!() not implemented for ", typeof(integrator))


function solutionstep(int::AbstractIntegrator, sol, args...)
    solstep = SolutionStep(problem(int), method(int))
    initialize!(solstep, sol, problem(int), args...)
    solstep
end

cache(::AbstractIntegrator, DT) = nothing
cache(::AbstractIntegrator) = nothing


"""
GeometricIntegrator

Collects all data structures needed by an integrator:

* `problem`: [`EquationProblem`](@ref) to solve
* `method`: integration method
* `cache`: temprary data structures needed by method
* `solver`: linear or nonlinear solver needed by method
* `iguess`: initial guess for implicit methods
* `projection`: optional projection method

Constructors:

```
GeometricIntegrator(problem::EquationProblem, method::GeometricMethod; solver = default_solver(method), iguess = default_iguess(method), projection = default_projection(method))
```

"""
struct GeometricIntegrator{
        MT <: GeometricMethod,
        PT <: AbstractProblem,
        CT <: CacheDict{PT,MT},
        ST <: Union{NonlinearSolver,SolverMethod},
        IT <: Extrapolation
    } <: AbstractIntegrator

    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
end

function GeometricIntegrator(
        problem::AbstractProblem,
        integratormethod::GeometricMethod,
        solvermethod::SolverMethod,
        iguess::Extrapolation;
        options = default_options(),
        method = initmethod(integratormethod, problem),
        caches = CacheDict(problem, method),
        solver = initsolver(solvermethod, options, method, caches)
    )
    GeometricIntegrator(problem, method, caches, solver, iguess)
end

function GeometricIntegrator(
        problem::AbstractProblem,
        method::GeometricMethod;
        solver = default_solver(method),
        initialguess = default_iguess(method),
        kwargs...
    )
    GeometricIntegrator(problem, method, solver, initialguess; kwargs...)
end

GeometricIntegrator(::AbstractProblem, ::Nothing, args...; kwargs...) = nothing

problem(int::GeometricIntegrator) = int.problem
method(int::GeometricIntegrator) = int.method
caches(int::GeometricIntegrator) = int.caches
solver(int::GeometricIntegrator) = int.solver
iguess(int::GeometricIntegrator) = int.iguess
initialguess(int::GeometricIntegrator) = int.iguess

cache(int::GeometricIntegrator, DT) = caches(int)[DT]
cache(int::GeometricIntegrator) = cache(int, datatype(problem(int)))
eachstage(int::GeometricIntegrator) = eachstage(method(int))
hasnullvector(int::GeometricIntegrator) = hasnullvector(method(int))
implicit_update(int::GeometricIntegrator) = implicit_update(method(int))
nconstraints(int::GeometricIntegrator) = nconstraints(problem(int))
Base.ndims(int::GeometricIntegrator) = ndims(problem(int))
nstages(int::GeometricIntegrator) = nstages(tableau(method(int)))
nlsolution(int::GeometricIntegrator) = nlsolution(cache(int))
nullvector(int::GeometricIntegrator) = nullvector(method(int))
tableau(int::GeometricIntegrator) = tableau(method(int))

equations(int::GeometricIntegrator) = functions(problem(int))
timestep(int::GeometricIntegrator) = timestep(problem(int))

initial_guess!(sol, history, params, ::GeometricIntegrator) = nothing

# Cache{DT}(int::Integrator) where {DT} = Cache{DT}(int.problem, int.method)

# cache(::GeometricIntegrator, DT) = nothing
# cache(::GeometricIntegrator) = nothing

# initialize!(::SolutionStep, ::AbstractProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
# initialize!(int::DeterministicIntegrator) = initialize!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))

# function residual!(b::AbstractVector, x::AbstractVector, int::DeterministicIntegrator)
#     residual!(b, x, current(solstep(int)), history(solstep(int)), parameters(solstep(int)), int)
# end

# components!(x::AbstractVector, int::DeterministicIntegrator) = components!(x, solstep(int), problem(int), method(int), caches(int))


# Parts of one integration step that are common to deterministic and stochastic equations.
# function integrate!(solstep::SolutionStep, problem::EquationProblem, method::GeometricMethod, caches::CacheDict, solver::Union{SolverMethod, NonlinearSolver}, iguess::Union{InitialGuess,Extrapolation})
function integrate!(solstep::SolutionStep, int::AbstractIntegrator)
    # reset solution step
    reset!(solstep, timestep(int))

    # compute initial guess
    initial_guess!(current(solstep), history(solstep), parameters(solstep), int)

    # integrate one initial condition for one time step
    integrate_step!(current(solstep), history(solstep), parameters(solstep), int)

    # copy internal variables from cache to solution step
    copy_internal_variables!(solstep, cache(int))

    # copy solver status to solution step
    # solver_status!(solver(int), solstep.internal[:solver])

    # take care of periodic solutions
    cut_periodic_solution!(solstep, periodicity(problem(int)))

    # update vector field for initial guess
    update_vector_fields!(solstep, problem(int))

    return solstep
end

"""
Solve for time steps n with n₁ ≤ n ≤ n₂.
```julia
integrate!(solution, integrator, n₁, n₂)
```
"""
function integrate!(sol::GeometricSolution, int::AbstractIntegrator, n₁::Int, n₂::Int)
    # check time steps range for consistency
    @assert n₁ ≥ 1
    @assert n₂ ≥ n₁
    @assert n₂ ≤ ntime(sol)

    # copy initial condition from solution to solutionstep and initialize
    solstep = solutionstep(int, sol[n₁-1])

    # loop over time steps
    for n in n₁:n₂
        # integrate one step and copy solution from cache to solution
        sol[n] = integrate!(solstep, int)

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
function integrate!(sol::GeometricSolution, int::AbstractIntegrator)
    integrate!(sol, int, 1, ntime(sol))
end


function integrate(integrator::AbstractIntegrator)
    solution = Solution(problem(integrator))
    integrate!(solution, integrator)
end

function integrate(problem::AbstractProblem, method::GeometricMethod; kwargs...)
    integrator = GeometricIntegrator(problem, method; kwargs...)
    integrate(integrator)
end

function integrate(problems::EnsembleProblem, method::GeometricMethod; kwargs...)
    solutions = Solution(problems)

    for (problem, solution) in zip(problems, solutions)
        integrator = GeometricIntegrator(problem, method; kwargs...)
        integrate!(solution, integrator)
    end

    return solutions
end





#*****************************************************************************#
# Legacy code (we need to see which parts of this are still needed)           #
#*****************************************************************************#

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
internal_variables(::Integrator) = NamedTuple()
```
Returns a `NamedTuple` containing all internal variables of an integrator that
shall be stored in an [`SolutionStep`](@ref). If there is no method for a
specific integrator implemented an empty `NamedTuple()` is returned.
"""
internal_variables(::GeometricIntegrator) = NamedTuple()
internal_variables(::Nothing) = NamedTuple()


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

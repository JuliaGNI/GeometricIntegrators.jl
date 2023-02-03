
struct NoSolver <: SolverMethod end
struct NoProjection <: ProjectionMethod end

default_solver(::GeometricMethod) = NoSolver()
default_iguess(::GeometricMethod) = NoInitialGuess()
default_projection(::GeometricMethod) = NoProjection()

initmethod(method::GeometricMethod) = method
initmethod(method::GeometricMethod, ::GeometricProblem) = initmethod(method)
initsolver(::SolverMethod, ::SolutionStep, ::GeometricProblem, ::GeometricMethod, ::CacheDict) = NoSolver()
# projection(::ProjectionMethod, ::GeometricProblem, ::GeometricMethod, ::CacheDict) = NoProjection()


"""
Integrator

Collects all data structures needed by an integrator:

* `problem`: [`GeometricProblem`](@ref) to solve
* `method`: integration method
* `cache`: temprary data structures needed by method
* `solver`: linear or nonlinear solver needed by method
* `iguess`: initial guess for implicit methods
* `projection`: optional projection method

Constructors:

```
Integrator(problem::GeometricProblem, method::GeometricMethod; solver = default_solver(method), iguess = default_iguess(method), projection = default_projection(method))
```

"""
struct Integrator{PT, MT, CT, ST, IT, PJT, SST}
    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
    projection::PJT
    solstep::SST

    function Integrator(problem::GeometricProblem, integratormethod::GeometricMethod, solvermethod::SolverMethod, iguess::Union{InitialGuess,Extrapolation}, prjctn::ProjectionMethod)
        method = initmethod(integratormethod, problem)
        caches = CacheDict(problem, method)
        solstep = SolutionStep(problem, method)
        solver = initsolver(solvermethod, solstep, problem, method, caches)
        # prjctn = projection(_projection, problem, method, caches)

        new{typeof(problem),
            typeof(method),
            typeof(caches),
            typeof(solver),
            typeof(iguess),
            typeof(prjctn),
            typeof(solstep)
           }(problem, method, caches, solver, iguess, prjctn, solstep)
    end
end

function Integrator(problem::GeometricProblem, method::GeometricMethod;
            solver = default_solver(method),
            initialguess = default_iguess(method),
            projection = default_projection(method))
    Integrator(problem, method, solver, initialguess, projection)
end

problem(int::Integrator) = int.problem
method(int::Integrator) = int.method
caches(int::Integrator) = int.caches
solver(int::Integrator) = int.solver
iguess(int::Integrator) = int.iguess
projct(int::Integrator) = int.projection
initialguess(int::Integrator) = int.iguess
projection(int::Integrator) = int.projection
solstep(int::Integrator) = int.solstep

cache(int::Integrator, DT) = caches(int)[DT]
cache(int::Integrator) = cache(int, datatype(solstep(int)))

GeometricBase.equations(int::Integrator) = functions(problem(int))
GeometricBase.timestep(int::Integrator) = timestep(problem(int))

# Cache{DT}(int::Integrator) where {DT} = Cache{DT}(int.problem, int.method)

initialize!(::SolutionStep, ::GeometricProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
initialize!(int::Integrator) = initialize!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))

initial_guess!(::SolutionStep, ::GeometricProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
initial_guess!(int::Integrator) = initial_guess!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))


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
function integrate(problem::GeometricProblem, method::GeometricMethod; kwargs...)
    integrator = Integrator(problem, method; kwargs...)
    solution = Solution(problem; kwargs...)
    integrate!(solution, integrator)
end


#*****************************************************************************#
# Integration functions for deterministic integrators                         #
#*****************************************************************************#

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
function integrate!(solstep::SolutionStep, problem::GeometricProblem, method::GeometricMethod, caches::CacheDict, solver::Union{SolverMethod, NonlinearSolver}, iguess::Union{InitialGuess,Extrapolation})
    # reset solution step
    reset!(solstep, timestep(problem))

    # compute initial guess
    initial_guess!(solstep, problem, method, caches, solver, iguess)

    # integrate one initial condition for one time step
    integrate_step!(solstep, problem, method, caches, solver)

    # take care of periodic solutions
    cut_periodic_solution!(solstep, periodicity(problem))

    return solstep
end

integrate!(int::Integrator) = integrate!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))

# Integrate equation for time steps n with n₁ ≤ n ≤ n₂.
function integrate!(sol::GeometricSolution, int::Integrator, n₁::Int, n₂::Int)
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
function integrate!(sol::GeometricSolution, int::Integrator)
    integrate!(sol, int, 1, ntime(sol))
end

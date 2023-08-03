"""
Integrator

Collects all data structures needed by an integrator:

* `problem`: [`EquationProblem`](@ref) to solve
* `method`: integration method
* `cache`: temprary data structures needed by method
* `solver`: linear or nonlinear solver needed by method
* `iguess`: initial guess for implicit methods
* `projection`: optional projection method

Constructors:

```
Integrator(problem::EquationProblem, method::GeometricMethod; solver = default_solver(method), iguess = default_iguess(method), projection = default_projection(method))
```

"""
struct Integrator{PT, MT, CT, ST, IT, SIT, SST} <: DeterministicIntegrator
    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
    subint::SIT
    solstep::SST

    function Integrator(
        problem::AbstractProblem, 
        integratormethod::GeometricMethod, 
        solvermethod::SolverMethod, 
        iguess::Union{InitialGuess,Extrapolation};
        method = initmethod(integratormethod, problem),
        caches = CacheDict(problem, method),
        subint = Integrator(problem, parent(method)),
        solstp = subint === nothing ? SolutionStep(problem, method) : solstep(subint),
        solver = initsolver(solvermethod, method, caches)
        )

        new{typeof(problem),
            typeof(method),
            typeof(caches),
            typeof(solver),
            typeof(iguess),
            typeof(subint),
            typeof(solstp)
           }(problem, method, caches, solver, iguess, subint, solstp)
    end
end

function Integrator(
    problem::AbstractProblem,
    method::GeometricMethod;
    solver = default_solver(method),
    initialguess = default_iguess(method), kwargs...)

    Integrator(problem, method, solver, initialguess; kwargs...)
end

Integrator(::AbstractProblem, ::Nothing, args...; kwargs...) = nothing

problem(int::Integrator) = int.problem
method(int::Integrator) = int.method
caches(int::Integrator) = int.caches
solver(int::Integrator) = int.solver
iguess(int::Integrator) = int.iguess
initialguess(int::Integrator) = int.iguess
subint(int::Integrator) = int.subint
solstep(int::Integrator) = int.solstep

cache(int::Integrator, DT) = caches(int)[DT]
cache(int::Integrator) = cache(int, datatype(solstep(int)))
eachstage(int::Integrator) = eachstage(method(int))
hasnullvector(int::Integrator) = hasnullvector(method(int))
implicit_update(int::Integrator) = implicit_update(method(int))
nconstraints(int::Integrator) = nconstraints(problem(int))
Base.ndims(int::Integrator) = ndims(problem(int))
nstages(int::Integrator) = nstages(tableau(method(int)))
nlsolution(int::Integrator) = nlsolution(cache(int))
nullvector(int::Integrator) = nullvector(method(int))
tableau(int::Integrator) = tableau(method(int))

GeometricBase.equations(int::Integrator) = functions(problem(int))
GeometricBase.timestep(int::Integrator) = timestep(problem(int))

# Cache{DT}(int::Integrator) where {DT} = Cache{DT}(int.problem, int.method)

initialize!(::SolutionStep, ::AbstractProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
initialize!(int::Integrator) = initialize!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))

initial_guess!(::SolutionStep, ::AbstractProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
initial_guess!(int::Integrator) = initial_guess!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))


function residual!(b::AbstractVector, x::AbstractVector, int::Integrator)
    residual!(b, x, solstep(int), problem(int), method(int), caches(int))
end

# components!(x, int) = components!(x, solstep(int), problem(int), method(int), caches(int))

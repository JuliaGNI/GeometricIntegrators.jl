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
struct GeometricIntegrator{PT, MT, CT, ST, IT, SIT, SST} <: DeterministicIntegrator
    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
    subint::SIT
    solstep::SST

    function GeometricIntegrator(
        problem::AbstractProblem, 
        integratormethod::GeometricMethod, 
        solvermethod::SolverMethod, 
        iguess::Union{InitialGuess,Extrapolation};
        method = initmethod(integratormethod, problem),
        caches = CacheDict(problem, method),
        subint = GeometricIntegrator(problem, parent(method)),
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

function GeometricIntegrator(
    problem::AbstractProblem,
    method::GeometricMethod;
    solver = default_solver(method),
    initialguess = default_iguess(method), kwargs...)

    GeometricIntegrator(problem, method, solver, initialguess; kwargs...)
end

GeometricIntegrator(::AbstractProblem, ::Nothing, args...; kwargs...) = nothing

problem(int::GeometricIntegrator) = int.problem
method(int::GeometricIntegrator) = int.method
caches(int::GeometricIntegrator) = int.caches
solver(int::GeometricIntegrator) = int.solver
iguess(int::GeometricIntegrator) = int.iguess
initialguess(int::GeometricIntegrator) = int.iguess
subint(int::GeometricIntegrator) = int.subint
solstep(int::GeometricIntegrator) = int.solstep

cache(int::GeometricIntegrator, DT) = caches(int)[DT]
cache(int::GeometricIntegrator) = cache(int, datatype(solstep(int)))
eachstage(int::GeometricIntegrator) = eachstage(method(int))
hasnullvector(int::GeometricIntegrator) = hasnullvector(method(int))
implicit_update(int::GeometricIntegrator) = implicit_update(method(int))
nconstraints(int::GeometricIntegrator) = nconstraints(problem(int))
Base.ndims(int::GeometricIntegrator) = ndims(problem(int))
nstages(int::GeometricIntegrator) = nstages(tableau(method(int)))
nlsolution(int::GeometricIntegrator) = nlsolution(cache(int))
nullvector(int::GeometricIntegrator) = nullvector(method(int))
tableau(int::GeometricIntegrator) = tableau(method(int))

GeometricBase.equations(int::GeometricIntegrator) = functions(problem(int))
GeometricBase.timestep(int::GeometricIntegrator) = timestep(problem(int))

# Cache{DT}(int::Integrator) where {DT} = Cache{DT}(int.problem, int.method)

initialize!(::SolutionStep, ::AbstractProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
initialize!(int::GeometricIntegrator) = initialize!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))

initial_guess!(::SolutionStep, ::AbstractProblem, ::GeometricMethod, ::CacheDict, ::Union{SolverMethod, NonlinearSolver}, ::Union{InitialGuess,Extrapolation}) = nothing
initial_guess!(int::GeometricIntegrator) = initial_guess!(solstep(int), problem(int), method(int), caches(int), solver(int), iguess(int))


function residual!(b::AbstractVector, x::AbstractVector, int::GeometricIntegrator)
    residual!(b, x, solstep(int), problem(int), method(int), caches(int))
end

# components!(x, int) = components!(x, solstep(int), problem(int), method(int), caches(int))

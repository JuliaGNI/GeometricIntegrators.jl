"""
A `ProjectionMethod` is an algorithm that is applied together with a geometric integrator
to enforce constraints which are not automatically satisfied by the integrator.
Examples include conservation of invariants like energy or the Dirac constraint in [`IODE`](@ref)s.
"""
abstract type ProjectionMethod <: GeometricMethod end

struct NoProjection <: ProjectionMethod end

projection(::GeometricMethod) = NoProjection()


"""
A `ProjectedMethod` consists of a [`ProjectionMethod`](@ref) and a [`GeometricMethod`](@ref).
"""
struct ProjectedMethod{PT <: ProjectionMethod, MT <: GeometricMethod} <: ProjectionMethod
    projection::PT
    method::MT

    function ProjectedMethod(proj::ProjectionMethod, method::GeometricMethod)
        _method = initmethod(method)
        new{typeof(proj), typeof(_method)}(proj, _method)
    end
end

projection(proj::ProjectedMethod) = proj.projection
Base.parent(proj::ProjectedMethod) = proj.method

initmethod(projection::ProjectedMethod) = projection
initmethod(projection::ProjectionMethod, method::GeometricMethod) = ProjectedMethod(projection, method)


@doc raw"""
The `ProjectionIntegrator` is the counterpart to the `GeometricIntegrator` for [`ProjectionMethod`](@ref)s.
"""
struct ProjectionIntegrator{
        MT <: ProjectionMethod,
        PT <: AbstractProblem,
        CT <: CacheDict{PT,MT},
        ST <: Union{NonlinearSolver,SolverMethod},
        IT <: Union{InitialGuess,Extrapolation},
        SIT <: AbstractIntegrator,
        SST <: SolutionStep
    } <: DeterministicIntegrator

    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
    subint::SIT
    solstep::SST
end

function ProjectionIntegrator(
        problem::AbstractProblem,
        projectionmethod::ProjectionMethod,
        solvermethod::SolverMethod,
        iguess::Union{InitialGuess,Extrapolation},
        subint::AbstractIntegrator;
        method = initmethod(projectionmethod, problem),
        caches = CacheDict(problem, method),
        solstp = solstep(subint),
        solver = initsolver(solvermethod, method, caches)
    )
    ProjectionIntegrator(problem, method, caches, solver, iguess, subint, solstp)
end

function ProjectionIntegrator(
        problem::AbstractProblem,
        projectionmethod::ProjectionMethod,
        solvermethod::SolverMethod,
        iguess::Union{InitialGuess,Extrapolation},
        parent_solvermethod::SolverMethod,
        parent_iguess::Union{InitialGuess,Extrapolation};
        kwargs...
    )
    subint = GeometricIntegrator(problem, parent(projectionmethod), parent_solvermethod, parent_iguess)
    ProjectionIntegrator(problem, projectionmethod, solvermethod, iguess, subint)
end

function GeometricIntegrator(
        problem::AbstractProblem,
        method::ProjectionMethod;
        solver = default_solver(method),
        initialguess = default_iguess(method),
        parent_solver = default_solver(parent(method)),
        parent_initialguess = default_iguess(parent(method)),
        kwargs...
    )
    ProjectionIntegrator(problem, method, solver, initialguess, parent_solver, parent_initialguess; kwargs...)
end


problem(int::ProjectionIntegrator) = int.problem
method(int::ProjectionIntegrator) = int.method
caches(int::ProjectionIntegrator) = int.caches
solver(int::ProjectionIntegrator) = int.solver
iguess(int::ProjectionIntegrator) = int.iguess
initialguess(int::ProjectionIntegrator) = int.iguess
subint(int::ProjectionIntegrator) = int.subint
solstep(int::ProjectionIntegrator) = int.solstep

cache(int::ProjectionIntegrator, DT) = caches(int)[DT]
cache(int::ProjectionIntegrator) = cache(int, datatype(solstep(int)))
nconstraints(int::ProjectionIntegrator) = nconstraints(problem(int))
nlsolution(int::ProjectionIntegrator) = nlsolution(cache(int))

Base.ndims(int::ProjectionIntegrator) = ndims(problem(int))

equations(int::ProjectionIntegrator) = functions(problem(int))
timestep(int::ProjectionIntegrator) = timestep(problem(int))

initialize!(int::ProjectionIntegrator) = initialize!(subint(int))
<<<<<<< HEAD
=======



function project!(solstep::SolutionStepODE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, timestep(problem))
end

function project!(solstep::SolutionStepDAE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, λ, timestep(problem))
end

function project!(solstep::SolutionStepPODE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, G, timestep(problem))
end

function project!(solstep::SolutionStepPDAE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, G, λ, timestep(problem))
end
>>>>>>> ed182cbc (Change order of type parameters of GeometricIntegrator and ProjectionIntegrator.)

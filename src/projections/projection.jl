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
struct ProjectedMethod{PT<:ProjectionMethod,MT<:GeometricMethod} <: ProjectionMethod
    projection::PT
    method::MT

    function ProjectedMethod(proj::ProjectionMethod, method::GeometricMethod)
        _method = initmethod(method)
        new{typeof(proj),typeof(_method)}(proj, _method)
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
    MT<:ProjectionMethod,
    PT<:AbstractProblem,
    CT<:CacheDict{PT,MT},
    ST<:Union{NonlinearSolver,SolverMethod},
    IT<:Extrapolation,
    SIT<:AbstractIntegrator
} <: AbstractIntegrator

    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
    subint::SIT
end

function ProjectionIntegrator(
    problem::AbstractProblem,
    projectionmethod::ProjectionMethod,
    solvermethod::SolverMethod,
    iguess::Extrapolation,
    subint::AbstractIntegrator;
    method=initmethod(projectionmethod, problem),
    caches=CacheDict(problem, method),
    options...
)
    solver = initsolver(solvermethod, method, caches; (length(options) == 0 ? default_options() : options)...)
    ProjectionIntegrator(problem, method, caches, solver, iguess, subint)
end

function ProjectionIntegrator(
    problem::AbstractProblem,
    projectionmethod::ProjectionMethod,
    solvermethod::SolverMethod,
    iguess::Extrapolation,
    parent_solvermethod::SolverMethod,
    parent_iguess::Extrapolation;
    kwargs...
)
    subint = GeometricIntegrator(problem, parent(projectionmethod), parent_solvermethod, parent_iguess; default_options(parent(projectionmethod))...)
    ProjectionIntegrator(problem, projectionmethod, solvermethod, iguess, subint; kwargs...)
end

function GeometricIntegrator(
    problem::AbstractProblem,
    method::ProjectionMethod;
    solver=default_solver(method),
    initialguess=default_iguess(method),
    parent_solver=default_solver(parent(method)),
    parent_initialguess=default_iguess(parent(method)),
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
# solstep(int::ProjectionIntegrator) = int.solstep

cache(int::ProjectionIntegrator, DT) = caches(int)[DT]
cache(int::ProjectionIntegrator) = cache(int, datatype(problem(int)))
nconstraints(int::ProjectionIntegrator) = nconstraints(problem(int))
nlsolution(int::ProjectionIntegrator) = nlsolution(cache(int))

Base.ndims(int::ProjectionIntegrator) = ndims(problem(int))

equations(int::ProjectionIntegrator) = functions(problem(int))
timestep(int::ProjectionIntegrator) = timestep(problem(int))

initialize!(int::ProjectionIntegrator) = initialize!(subint(int))

# const DAEProjectionIntegrator{MT} = ProjectionIntegrator{MT, <:DAEProblem} where {{MT <: ProjectionMethod}}
# const IODEProjectionIntegrator{MT} = ProjectionIntegrator{MT, <:IODEProblem} where {{MT <: ProjectionMethod}}
# const LODEProjectionIntegrator{MT} = ProjectionIntegrator{MT, <:LODEProblem} where {{MT <: ProjectionMethod}}


function project!(sol, U::AbstractVector, G::AbstractVector, int::ProjectionIntegrator{<:ProjectionMethod,<:DAEProblem})
    sol.q .+= timestep(int) .* U
end

function project!(sol, U::AbstractVector, G::AbstractVector, int::ProjectionIntegrator{<:ProjectionMethod,<:Union{IODEProblem,LODEProblem,PDAEProblem,HDAEProblem}})
    sol.q .+= timestep(int) .* U
    sol.p .+= timestep(int) .* G
end

# function project!(sol, problem::EquationProblem, method::ProjectionMethod, cache::ProjectionCache)
#     project(solstep, problem, method, cache.U, cache.G, cache.λ)
# end

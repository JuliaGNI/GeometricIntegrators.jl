"""
    linear_solver_defaults(solvermethod)

The keyword arguments that set the linear solver of a nonlinear solver built by this package.
A caller's keyword arguments are splatted after them, so a caller who names
`linear_solver_method` overrides the default.

The linear solver of `Newton`, `QuasiNewton` and `DogLeg` is SimpleSolvers' pure-Julia `LU()`.
From SimpleSolvers 0.13 the default is `LapackLU`, whose different order of operations turns a
numerically singular stage system from a `SingularException` into a silently wrong solution.
`Picard` takes no linear solver, so it gets no keyword.
"""
linear_solver_defaults(::SolverMethod) = NamedTuple()
function linear_solver_defaults(::Union{Newton, DogLeg})
    (linear_solver_method = SimpleSolvers.LU(),)
end

# every method type of this package whose supertype belongs to another package
const OwnMethod = Union{
    AbstractSplittingMethod, AbstractCompositionMethod, Composition, Splitting, ExactSolution,
    RKMethod, PRKMethod, PGLRK, FLRK,
    VIMethod, DVIMethod, CGVIMethod, DGVIMethod, HPIMethod, DiscreteEulerLagrange,
    ProjectedMethod, StandardProjection, MidpointProjection, SymmetricProjection,
    VariationalProjection, InternalStageProjection, LegendreProjection, SecondaryProjection
}

function initsolver(
        solvermethod::NonlinearSolverMethod, method::OwnMethod, caches::CacheDict;
        kwargs...)
    invoke(initsolver, Tuple{NonlinearSolverMethod, GeometricMethod, CacheDict},
        solvermethod, method, caches; linear_solver_defaults(solvermethod)..., kwargs...)
end

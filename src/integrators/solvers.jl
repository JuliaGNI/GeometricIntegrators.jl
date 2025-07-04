
struct NoSolver <: SolverMethod end

default_linesearch(method=nothing) = Backtracking()

default_options(method=nothing) = (
    min_iterations=1,
    f_abstol=8eps(),
    linesearch=default_linesearch(method),
    # verbosity=2,
)

initsolver(::SolverMethod, ::GeometricMethod, ::CacheDict; kwargs...) = NoSolver()

# create nonlinear solver
function initsolver(::NewtonMethod, method::GeometricMethod, caches::CacheDict; kwargs...)
    x = zero(nlsolution(caches))
    y = zero(nlsolution(caches))
    NewtonSolver(x, residual!, y; kwargs...)
end

# residual!(y, x, parameters::Union{Tuple,NamedTuple}) = residual!(y, x, parameters...)
residual!(y, x, parameters) = residual!(y, x, parameters...) # TODO: This is a workaround for a SimpleSolvers limitation. Should be removed.

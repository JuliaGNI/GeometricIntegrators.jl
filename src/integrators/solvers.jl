
struct NoSolver <: SolverMethod end

default_options() = Options(
    min_iterations = 1,
    x_abstol = 8eps(),
    f_abstol = 8eps(),
)

initsolver(::SolverMethod, ::Options, ::GeometricMethod, ::CacheDict) = NoSolver()

# create nonlinear solver
function initsolver(::NewtonMethod, config::Options, ::GeometricMethod, caches::CacheDict; kwargs...)
    x = zero(nlsolution(caches))
    y = zero(nlsolution(caches))
    NewtonSolver(x, y; linesearch = Backtracking(), config = config, kwargs...)
end

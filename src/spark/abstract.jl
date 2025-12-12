
abstract type HSPARKMethod <: HDAEMethod end
abstract type ISPARKMethod <: IDAEMethod end
abstract type LSPARKMethod <: LDAEMethod end
abstract type PSPARKMethod <: PDAEMethod end

const AbstractSPARKMethod = Union{HSPARKMethod,ISPARKMethod,LSPARKMethod,PSPARKMethod}
const AbstractSPARKProblem{DT<:Number,TT<:Real} =
    Union{HDAEProblem{DT,TT},
        IDAEProblem{DT,TT},
        LDAEProblem{DT,TT},
        PDAEProblem{DT,TT}}

GeometricIntegratorsBase.default_iguess(::AbstractSPARKMethod) = HermiteExtrapolation()
GeometricIntegratorsBase.default_solver(::AbstractSPARKMethod) = Newton()
GeometricIntegratorsBase.default_options(::AbstractSPARKMethod) = (
    min_iterations=1,
    x_suctol=2eps(),
    f_abstol=8eps(),
    f_suctol=2eps(),
)


nstages(method::AbstractSPARKMethod) = nstages(tableau(method))
pstages(method::AbstractSPARKMethod) = pstages(tableau(method))
eachstage(method::AbstractSPARKMethod) = eachstage(tableau(method))
hasnullvector(method::AbstractSPARKMethod) = hasnullvector(tableau(method))

# create nonlinear solver
function initsolver(::Newton, method::AbstractSPARKMethod, caches::CacheDict; kwargs...)
    NewtonSolver(zero(nlsolution(caches)), residual!, zero(nlsolution(caches)); linesearch=Backtracking(), kwargs...)
end

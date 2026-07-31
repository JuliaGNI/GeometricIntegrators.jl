
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
# The SPARK residual settles at a round-off floor between ~4e-16 (well-conditioned tableaux
# such as SLRK and VSPARK(SPARKLobABC)) and ~3e-15 (the marginal VSPARK(SPARKLobABD(4))).
# An absolute tolerance one order of magnitude above that floor lets these solves converge
# silently, rather than stalling against the unreachable framework default f_abstol = 0 and
# emitting a stagnation warning at the finest steps. 8e-15 was measured to leave the computed
# solutions unchanged while clearing the SPARKLobABD(4) warning; see docs/src/audit.md (fourth
# pass). x_suctol/f_suctol are left at the SimpleSolvers default (2eps) and so need no restatement.
GeometricIntegratorsBase.default_options(::AbstractSPARKMethod) = (
    min_iterations = 1,
    f_abstol = 8e-15,
)


nstages(method::AbstractSPARKMethod) = nstages(tableau(method))
pstages(method::AbstractSPARKMethod) = pstages(tableau(method))
eachstage(method::AbstractSPARKMethod) = eachstage(tableau(method))
hasnullvector(method::AbstractSPARKMethod) = hasnullvector(tableau(method))

"""
    nullvectorsize(problem, method)

Number of unknowns the null-vector multiplier `μ` contributes to the nonlinear system:
one component per degree of freedom when the tableau carries a null vector `d`, zero
otherwise.

Every `solversize` method in this submodule adds this term, so `solversize` is the full
length of the solver's unknown vector and the cache constructor needs no correction of
its own. Keep it that way — the two used to be computed in different files, and nothing
checked that they agreed (see the `SPARK solver size` testset).
"""
nullvectorsize(problem::AbstractSPARKProblem, method::AbstractSPARKMethod) =
    hasnullvector(method) ? length(vec(initial_conditions(problem).q)) : 0

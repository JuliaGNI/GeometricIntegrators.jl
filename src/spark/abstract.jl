
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
# No `default_options` override: SPARK relies on the GeometricIntegratorsBase default
# `f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))`. SPARK's stage
# systems are high-dimensional (solversize ≫ 8), so `solversize · eps` sits above their
# round-off floor (~4e-16 for well-conditioned tableaux, ~3e-15 for the marginal
# VSPARK(SPARKLobABD(4))), which is why the flat 8e-15 override is no longer needed.


nstages(method::AbstractSPARKMethod) = nstages(tableau(method))
pstages(method::AbstractSPARKMethod) = pstages(tableau(method))
eachstage(method::AbstractSPARKMethod) = eachstage(tableau(method))
hasnullvector(method::AbstractSPARKMethod) = hasnullvector(tableau(method))

"""
    nullvectorsize(method, problem)

Number of unknowns the null-vector multiplier `μ` contributes to the nonlinear system:
one component per degree of freedom when the tableau carries a null vector `d`, zero
otherwise.

Every `solversize` method in this submodule adds this term, so `solversize` is the full
length of the solver's unknown vector and the cache constructor needs no correction of
its own. Keep it that way — the two used to be computed in different files, and nothing
checked that they agreed (see the `SPARK solver size` testset).
"""
nullvectorsize(method::AbstractSPARKMethod, problem::AbstractSPARKProblem) =
    hasnullvector(method) ? length(vec(initial_conditions(problem).q)) : 0

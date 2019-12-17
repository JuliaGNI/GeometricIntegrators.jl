
"""
Abstract atomistic or single-step solution.
"""
abstract type AtomisticSolution{dType <: Number, tType <: Real} end


"Create AtomisticSolution for ODE."
function AtomisticSolution(equation::AbstractEquationODE{DT,TT}) where {DT,TT}
    AtomisticSolutionODE{DT,TT}(ndims(equation))
end

"Create AtomisticSolution for partitioned ODE."
function AtomisticSolution(equation::AbstractEquationPODE{DT,TT}) where {DT,TT}
    AtomisticSolutionPODE{DT,TT}(ndims(equation))
end

"Create AtomisticSolution for DAE."
function AtomisticSolution(equation::AbstractEquationDAE{DT,TT}) where {DT,TT}
    AtomisticSolutionDAE{DT,TT}(ndims(equation))
end

"Create AtomisticSolution for partitioned DAE."
function AtomisticSolution(equation::AbstractEquationPDAE{DT,TT}) where {DT,TT}
    AtomisticSolutionPDAE{DT,TT}(ndims(equation))
end

"Print error for AtomisticSolutions of equations not implemented, yet."
function AtomisticSolution(equation::Equation)
    error("No AtomisticSolution found for equation ", equation)
end

"Copy solution from atomistic solution to solution object."
function copy_solution!(sol::Solution, asol::AtomisticSolution, n, m)
    copy_solution!(sol, get_solution(asol)..., n, m)
end

function CommonFunctions.cut_periodic_solution!(asol::AtomisticSolution{DT}, periodicity::Vector{DT}) where {DT}
    @assert length(asol.q) == length(periodicity)
    for k in eachindex(asol.q, asol.q̅, periodicity)
        if periodicity[k] ≠ 0
            while asol.q[k] < 0
                asol.q[k], asol.q̃[k] = compensated_summation(+periodicity[k], asol.q[k], asol.q̃[k])
                asol.q̅[k] += periodicity[k]
            end
            while asol.q[k] ≥ periodicity[k]
                asol.q[k], asol.q̃[k] = compensated_summation(-periodicity[k], asol.q[k], asol.q̃[k])
                asol.q̅[k] -= periodicity[k]
            end
        end
    end
end


"""
Abstract atomistic or single-step solution.
"""
abstract type AtomicSolution{dType <: Number, tType <: Real} end


"Create AtomicSolution for ODE."
function AtomicSolution(equation::AbstractEquationODE{DT,TT}) where {DT,TT}
    AtomicSolutionODE(DT, TT, ndims(equation))
end

"Create AtomicSolution for partitioned ODE."
function AtomicSolution(equation::AbstractEquationPODE{DT,TT}) where {DT,TT}
    AtomicSolutionPODE(DT, TT, ndims(equation))
end

"Create AtomicSolution for DAE."
function AtomicSolution(equation::AbstractEquationDAE{DT,TT}) where {DT,TT}
    AtomicSolutionDAE(DT, TT, ndims(equation), equation.m)
end

"Create AtomicSolution for partitioned DAE."
function AtomicSolution(equation::AbstractEquationPDAE{DT,TT}) where {DT,TT}
    AtomicSolutionPDAE(DT, TT, ndims(equation), equation.m)
end

"Create AtomicSolution for SDE."
function AtomicSolution(equation::AbstractEquationSDE{DT,TT}) where {DT,TT}
    AtomicSolutionSDE(DT, TT, ndims(equation), equation.m)
end

"Create AtomicSolution for PSDE."
function AtomicSolution(equation::AbstractEquationPSDE{DT,TT}) where {DT,TT}
    AtomicSolutionPSDE(DT, TT, ndims(equation), equation.m)
end

"Print error for AtomicSolutions of equations not implemented, yet."
function AtomicSolution(equation::Equation)
    error("No AtomicSolution found for equation ", equation)
end

"Copy solution from atomistic solution to solution object."
function copy_solution!(sol::Solution, asol::AtomicSolution, n, m)
    copy_solution!(sol, get_solution(asol)..., n, m)
end

function CommonFunctions.cut_periodic_solution!(asol::AtomicSolution{DT}, periodicity::Vector{DT}) where {DT}
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

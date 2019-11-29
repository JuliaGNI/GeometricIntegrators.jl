
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

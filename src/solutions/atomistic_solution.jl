
"""
Abstract atomistic or single-step solution.
"""
abstract type AtomisticSolution{dType <: Number, tType <: Real} end


"Create AtomisticSolution for ODE."
function AtomisticSolution(equation::AbstractEquationODE)
    AtomisticSolutionODE()
end

"Create AtomisticSolution for partitioned ODE."
function AtomisticSolution(equation::AbstractEquationPODE)
    AtomisticSolutionPODE()
end

"Create AtomisticSolution for DAE."
function AtomisticSolution(equation::AbstractEquationDAE)
    AtomisticSolutionDAE()
end

"Create AtomisticSolution for partitioned DAE."
function AtomisticSolution(equation::AbstractEquationPDAE)
    AtomisticSolutionPDAE()
end

"Print error for AtomisticSolutions of equations not implemented, yet."
function AtomisticSolution(equation::Equation)
    error("No AtomisticSolution found for equation ", equation)
end

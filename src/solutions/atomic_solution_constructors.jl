
"Create AtomicSolution for ODE."
function AtomicSolution(equation::AbstractEquationODE)
    AtomicSolutionODE(equation.t₀, equation.q₀[begin])
end

"Create AtomicSolution for ODE."
function AtomicSolution(solution::SolutionODE)
    AtomicSolutionODE(get_initial_conditions(solution, 1)...)
end

"Create AtomicSolution for partitioned ODE."
function AtomicSolution(equation::AbstractEquationPODE)
    AtomicSolutionPODE(equation.t₀, equation.q₀[begin], equation.p₀[begin])
end

"Create AtomicSolution for partitioned ODE."
function AtomicSolution(solution::SolutionPODE)
    AtomicSolutionPODE(get_initial_conditions(solution, 1)...)
end

"Create AtomicSolution for DAE."
function AtomicSolution(equation::AbstractEquationDAE)
    AtomicSolutionDAE(equation.t₀, equation.q₀[begin], equation.λ₀[begin])
end

"Create AtomicSolution for DAE."
function AtomicSolution(solution::SolutionDAE)
    AtomicSolutionDAE(get_initial_conditions(solution, 1)...)
end

"Create AtomicSolution for partitioned DAE."
function AtomicSolution(equation::AbstractEquationPDAE)
    AtomicSolutionPDAE(equation.t₀, equation.q₀[begin], equation.p₀[begin], equation.λ₀[begin])
end

"Create AtomicSolution for partitioned DAE."
function AtomicSolution(solution::SolutionPDAE)
    AtomicSolutionPDAE(get_initial_conditions(solution, 1)...)
end

"Print error for AtomicSolutions of equations not implemented, yet."
function AtomicSolution(equation::Equation)
    error("No AtomicSolution found for equation ", equation)
end

"Print error for AtomicSolutions of solution not implemented, yet."
function AtomicSolution(solution::Solution)
    error("No AtomicSolution found for solution ", solution)
end

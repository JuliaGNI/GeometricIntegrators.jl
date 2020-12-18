
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

"Create AtomicSolution for SDE."
function AtomicSolution(equation::AbstractEquationSDE{DT,TT}) where {DT,TT}
    AtomicSolutionSDE(equation.t₀, equation.q₀[begin], zeros(DT,equation.m), zeros(DT,equation.m))
end

"Create AtomicSolution for SDE."
function AtomicSolution(solution::SolutionSDE{AT,TT}) where {DT, TT, AT <: AbstractArray{DT}}
    AtomicSolutionSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm))
end

"Create AtomicSolution for partitioned SDE."
function AtomicSolution(equation::AbstractEquationPSDE{DT,TT}) where {DT,TT}
    AtomicSolutionPSDE(equation.t₀, equation.q₀[begin], equation.p₀[begin], zeros(DT,equation.m), zeros(DT,equation.m))
end

"Create AtomicSolution for partitioned SDE."
function AtomicSolution(solution::SolutionPSDE{AT,TT}) where {DT, TT, AT <: AbstractArray{DT}}
    AtomicSolutionPSDE(get_initial_conditions(solution, 1)..., zeros(DT,solution.nm), zeros(DT,solution.nm))
end

"Print error for AtomicSolutions of equations not implemented, yet."
function AtomicSolution(equation::Equation)
    error("No AtomicSolution found for equation ", equation)
end

"Print error for AtomicSolutions of solution not implemented, yet."
function AtomicSolution(solution::Solution)
    error("No AtomicSolution found for solution ", solution)
end

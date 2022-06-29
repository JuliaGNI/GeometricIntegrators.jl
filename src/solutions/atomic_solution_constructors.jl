
# Create AtomicSolution for ODE.
function AtomicSolution(problem::AbstractProblemODE)
    AtomicSolutionODE(problem.tspan[begin], problem.ics.q)
end

# Create AtomicSolution for ODE.
function AtomicSolution(solution::SolutionODE)
    AtomicSolutionODE(get_initial_conditions(solution, 1)...)
end

# Create AtomicSolution for partitioned ODE.
function AtomicSolution(problem::AbstractProblemPODE)
    AtomicSolutionPODE(problem.tspan[begin], problem.ics.q, problem.ics.p)
end

# Create AtomicSolution for partitioned ODE.
function AtomicSolution(solution::SolutionPODE)
    AtomicSolutionPODE(get_initial_conditions(solution, 1)...)
end

# Create AtomicSolution for DAE.
function AtomicSolution(problem::AbstractProblemDAE)
    AtomicSolutionDAE(problem.tspan[begin], problem.ics.q, problem.ics.λ)
end

# Create AtomicSolution for DAE.
function AtomicSolution(solution::SolutionDAE)
    AtomicSolutionDAE(get_initial_conditions(solution, 1)...)
end

# Create AtomicSolution for partitioned DAE.
function AtomicSolution(problem::AbstractProblemPDAE)
    AtomicSolutionPDAE(problem.tspan[begin], problem.ics.q, problem.ics.p, problem.ics.λ)
end

# Create AtomicSolution for partitioned DAE.
function AtomicSolution(solution::SolutionPDAE)
    AtomicSolutionPDAE(get_initial_conditions(solution, 1)...)
end

# Print error for AtomicSolutions of problem types not implemented, yet.
function AtomicSolution(problem::GeometricProblem)
    error("No AtomicSolution found for problem type ", typeof(problem))
end

# Print error for AtomicSolutions of solution not implemented, yet.
function AtomicSolution(solution::Solution)
    error("No AtomicSolution found for solution ", solution)
end

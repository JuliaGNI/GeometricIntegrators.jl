
# Create AtomicSolution for ODE.
function AtomicSolution(problem::AbstractProblemODE, internal::NamedTuple = NamedTuple())
    AtomicSolutionODE(initial_conditions(problem)..., internal)
end

# Create AtomicSolution for ODE.
function AtomicSolution(solution::SolutionODE, internal::NamedTuple = NamedTuple())
    AtomicSolutionODE(solution[0]..., internal)
end

# Create AtomicSolution for partitioned ODE.
function AtomicSolution(problem::AbstractProblemPODE, internal::NamedTuple = NamedTuple())
    ics = initial_conditions(problem)
    AtomicSolutionPODE(ics.t, ics.q, ics.p, internal)
end

# Create AtomicSolution for partitioned ODE.
function AtomicSolution(solution::SolutionPODE, internal::NamedTuple = NamedTuple())
    AtomicSolutionPODE(solution[0].t, solution[0].q, solution[0].p, internal)
end

# Create AtomicSolution for DAE.
function AtomicSolution(problem::AbstractProblemDAE, internal::NamedTuple = NamedTuple())
    AtomicSolutionDAE(initial_conditions(problem)..., internal)
end

# Create AtomicSolution for DAE.
function AtomicSolution(solution::SolutionDAE, internal::NamedTuple = NamedTuple())
    AtomicSolutionDAE(solution[0]..., internal)
end

# Create AtomicSolution for partitioned DAE.
function AtomicSolution(problem::AbstractProblemPDAE, internal::NamedTuple = NamedTuple())
    ics = initial_conditions(problem)
    AtomicSolutionPDAE(ics.t, ics.q, ics.p, ics.λ, internal)
end

# Create AtomicSolution for partitioned DAE.
function AtomicSolution(solution::SolutionPDAE, internal::NamedTuple = NamedTuple())
    AtomicSolutionPDAE(solution[0].t, solution[0].q, solution[0].p, solution[0].λ, internal)
end

# Print error for AtomicSolutions of problem types not implemented, yet.
function AtomicSolution(problem::GeometricProblem, args...)
    error("No AtomicSolution found for problem type ", typeof(problem))
end

# Print error for AtomicSolutions of solution not implemented, yet.
function AtomicSolution(solution::AbstractSolution, args...)
    error("No AtomicSolution found for solution ", solution)
end

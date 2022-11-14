
# Create SolutionStep for an ODEProblem.
function SolutionStep(problem::AbstractProblemODE, internal::NamedTuple = NamedTuple())
    SolutionStepODE(initial_conditions(problem)..., internal)
end

# Create SolutionStep for an ODEProblem.
function SolutionStep(solution::SolutionODE, internal::NamedTuple = NamedTuple())
    SolutionStepODE(solution[0]..., internal)
end

# Create SolutionStep for a PODEProblem.
function SolutionStep(problem::AbstractProblemPODE, internal::NamedTuple = NamedTuple())
    ics = initial_conditions(problem)
    SolutionStepPODE(ics.t, ics.q, ics.p, internal)
end

# Create SolutionStep for a PODEProblem.
function SolutionStep(solution::SolutionPODE, internal::NamedTuple = NamedTuple())
    SolutionStepPODE(solution[0].t, solution[0].q, solution[0].p, internal)
end

# Create SolutionStep for a DAEProblem.
function SolutionStep(problem::AbstractProblemDAE, internal::NamedTuple = NamedTuple())
    SolutionStepDAE(initial_conditions(problem)..., internal)
end

# Create SolutionStep for a DAEProblem.
function SolutionStep(solution::SolutionDAE, internal::NamedTuple = NamedTuple())
    SolutionStepDAE(solution[0]..., internal)
end

# Create SolutionStep for a PDAEProblem.
function SolutionStep(problem::AbstractProblemPDAE, internal::NamedTuple = NamedTuple())
    ics = initial_conditions(problem)
    SolutionStepPDAE(ics.t, ics.q, ics.p, ics.λ, internal)
end

# Create SolutionStep for a PDAEProblem.
function SolutionStep(solution::SolutionPDAE, internal::NamedTuple = NamedTuple())
    SolutionStepPDAE(solution[0].t, solution[0].q, solution[0].p, solution[0].λ, internal)
end

# Print error for SolutionSteps of problem types not implemented, yet.
function SolutionStep(problem::GeometricProblem, args...)
    error("No SolutionStep found for problem type ", typeof(problem))
end

# Print error for SolutionSteps of solution not implemented, yet.
function SolutionStep(solution::AbstractSolution, args...)
    error("No SolutionStep found for solution ", solution)
end

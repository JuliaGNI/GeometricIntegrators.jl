
# Create SolutionStep for an ODEProblem.
function SolutionStep(problem::AbstractProblemODE; kwargs...)
    SolutionStepODE(initial_conditions(problem)...; kwargs...)
end

# Create SolutionStep for an ODEProblem.
function SolutionStep(solution::SolutionODE; kwargs...)
    SolutionStepODE(solution[0]...; kwargs...)
end

# Create SolutionStep for a PODEProblem.
function SolutionStep(problem::AbstractProblemPODE; kwargs...)
    ics = initial_conditions(problem)
    SolutionStepPODE(ics.t, ics.q, ics.p; kwargs...)
end

# Create SolutionStep for a PODEProblem.
function SolutionStep(solution::SolutionPODE; kwargs...)
    SolutionStepPODE(solution[0].t, solution[0].q, solution[0].p; kwargs...)
end

# Create SolutionStep for a DAEProblem.
function SolutionStep(problem::AbstractProblemDAE; kwargs...)
    SolutionStepDAE(initial_conditions(problem)...; kwargs...)
end

# Create SolutionStep for a DAEProblem.
function SolutionStep(solution::SolutionDAE; kwargs...)
    SolutionStepDAE(solution[0]...; kwargs...)
end

# Create SolutionStep for a PDAEProblem.
function SolutionStep(problem::AbstractProblemPDAE; kwargs...)
    ics = initial_conditions(problem)
    SolutionStepPDAE(ics.t, ics.q, ics.p, ics.λ; kwargs...)
end

# Create SolutionStep for a PDAEProblem.
function SolutionStep(solution::SolutionPDAE; kwargs...)
    SolutionStepPDAE(solution[0].t, solution[0].q, solution[0].p, solution[0].λ; kwargs...)
end

function SolutionStep(problem::GeometricProblem, method::GeometricMethod)
    SolutionStep(problem, internal(method))
end

function SolutionStep(solution::AbstractSolution, method::GeometricMethod)
    SolutionStep(solution, internal(method))
end

# Print error for SolutionSteps of problem types not implemented, yet.
function SolutionStep(problem::GeometricProblem, args...; kwargs...)
    error("No SolutionStep found for problem type ", typeof(problem))
end

# Print error for SolutionSteps of solution not implemented, yet.
function SolutionStep(solution::AbstractSolution, args...; kwargs...)
    error("No SolutionStep found for solution ", solution)
end


# Create SolutionStep for an ODEProblem.
function SolutionStep(problem::Union{ODEProblem, SODEProblem, SubstepProblem}, extrap::Extrapolation = default_extrapolation(); kwargs...)
    solstep = SolutionStepODE(initial_conditions(problem)...; kwargs...)
    initialize!(solstep, problem, extrap)
    return solstep
end

# Create SolutionStep for an ODEProblem.
function SolutionStep(solution::SolutionODE; kwargs...)
    SolutionStepODE(solution[0]...; kwargs...)
end

# Create SolutionStep for a PODEProblem.
function SolutionStep(problem::Union{PODEProblem, HODEProblem, IODEProblem, LODEProblem}, extrap::Extrapolation = default_extrapolation(); kwargs...)
    ics = initial_conditions(problem)
    solstep = SolutionStepPODE(ics.t, ics.q, ics.p; kwargs...)
    initialize!(solstep, problem, extrap)
    return solstep
end

# Create SolutionStep for a PODEProblem.
function SolutionStep(solution::SolutionPODE; kwargs...)
    SolutionStepPODE(solution[0].t, solution[0].q, solution[0].p; kwargs...)
end

# Create SolutionStep for a DAEProblem.
function SolutionStep(problem::DAEProblem, extrap::Extrapolation = default_extrapolation(); kwargs...)
    solstep = SolutionStepDAE(initial_conditions(problem)...; kwargs...)
    initialize!(solstep, problem, extrap)
    return solstep
end

# Create SolutionStep for a DAEProblem.
function SolutionStep(solution::SolutionDAE; kwargs...)
    SolutionStepDAE(solution[0]...; kwargs...)
end

# Create SolutionStep for a PDAEProblem.
function SolutionStep(problem::Union{PDAEProblem, HDAEProblem, IDAEProblem, LDAEProblem}, extrap::Extrapolation = default_extrapolation(); kwargs...)
    ics = initial_conditions(problem)
    solstep = SolutionStepPDAE(ics.t, ics.q, ics.p, ics.λ; kwargs...)
    initialize!(solstep, problem, extrap)
    return solstep
end

# Create SolutionStep for a PDAEProblem.
function SolutionStep(solution::SolutionPDAE; kwargs...)
    SolutionStepPDAE(solution[0].t, solution[0].q, solution[0].p, solution[0].λ; kwargs...)
end

function SolutionStep(problem::AbstractProblem, method::GeometricMethod, extrap::Extrapolation = default_extrapolation())
    SolutionStep(problem, extrap; internal=internal_variables(method))
end

function SolutionStep(solution::AbstractSolution, method::GeometricMethod)
    SolutionStep(solution; internal=internal_variables(method))
end

# Print error for SolutionSteps of problem types not implemented, yet.
function SolutionStep(problem::AbstractProblem, args...; kwargs...)
    error("No SolutionStep found for problem type ", typeof(problem))
end

# Print error for SolutionSteps of solution not implemented, yet.
function SolutionStep(solution::AbstractSolution, args...; kwargs...)
    error("No SolutionStep found for solution ", solution)
end

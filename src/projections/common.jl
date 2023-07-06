
function project!(solstep::SolutionStepODE, problem::GeometricProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, timestep(problem))
end

function project!(solstep::SolutionStepDAE, problem::GeometricProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, λ, timestep(problem))
end

function project!(solstep::SolutionStepPODE, problem::GeometricProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, G, timestep(problem))
end

function project!(solstep::SolutionStepPDAE, problem::GeometricProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, G, λ, timestep(problem))
end

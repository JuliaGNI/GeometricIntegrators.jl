
function project!(solstep::SolutionStepODE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, timestep(problem))
end

function project!(solstep::SolutionStepDAE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, λ, timestep(problem))
end

function project!(solstep::SolutionStepPODE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, G, timestep(problem))
end

function project!(solstep::SolutionStepPDAE, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    Solutions.update!(solstep, U, G, λ, timestep(problem))
end


function _projection_weights(RU, RG, R∞=1)
    DT = promote_type(eltype(RU), eltype(RG), typeof(R∞))
    RU = Vector{DT}([RU[1], R∞ * RU[2]])
    RG = Vector{DT}([RG[1], R∞ * RG[2]])
    return (DT, RU, RG)
end


function project!(solstep::SolutionStep{ODE}, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, timestep(problem))
end

function project!(solstep::SolutionStep{DAE}, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, λ, timestep(problem))
end

function project!(solstep::SolutionStep{PODE}, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, G, timestep(problem))
end

function project!(solstep::SolutionStep{PDAE}, problem::EquationProblem, ::ProjectionMethod, U, G, λ)
    update!(solstep, U, G, λ, timestep(problem))
end

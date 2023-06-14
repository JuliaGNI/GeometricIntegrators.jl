
function initialguess!(t̄, q̄, t, q, q̇, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, t, q, problem, extrap)
    functions(problem).v(q̇, t, q)
    (t = t, q = q, v = q̇)
end

function initialguess!(t, q, q̇, solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    initialguess!(previous(solstep)..., t, q, q̇, problem, extrap)
end

function initialguess!(t̄, q̄, t, q, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, t, q, problem, extrap)
    (t = t, q = q)
end

function initialguess!(t, q, solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    initialguess!(previous(solstep)..., t, q, problem, extrap)
end



function initialguess!(t̄, q̄, p̄, t, q, p, q̇, ṗ, problem::Union{PODEProblem,HODEProblem}, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, p̄, t, q, p, problem, extrap)
    functions(problem).v(q̇, t, q, p)
    functions(problem).f(ṗ, t, q, p)
    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t̄, q̄, p̄, t, q, p, q̇, ṗ, problem::Union{IODEProblem,LODEProblem}, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, p̄, t, q, p, problem, extrap)
    functions(problem).v̄(q̇, t, q)
    functions(problem).f̄(ṗ, t, q, q̇)
    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepPODE, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    initialguess!(solstep.t̄[1], solstep.q̄[1], solstep.p̄[1], t, q, p, q̇, ṗ, problem, extrap)
end

function initialguess!(t̄, q̄, p̄, t, q, p, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    # extrapolate!(t̄, q̄, p̄, t, q, p, problem, extrap)
    extrapolate!(t̄, q̄, t, q, problem, extrap)
    (t = t, q = q, p = p)
end

function initialguess!(t, q, p, solstep::SolutionStepPODE, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    initialguess!(solstep.t̄[1], solstep.q̄[1], solstep.p̄[1], t, q, p, problem, extrap)
end

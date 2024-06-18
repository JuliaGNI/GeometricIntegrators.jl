
function initialguess!(t̄, q̄, t, q, q̇, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    sol = (t = t, q = q)
    hist = (t = [t̄], q = [q̄])
    extrapolate!(sol, hist, problem, extrap)
    functions(problem).v(q̇, t, q, parameters(problem))
    (t = t, q = q, v = q̇)
end

function initialguess!(t, q, q̇, solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    initialguess!(previous(solstep)..., t, q, q̇, problem, extrap)
end

function initialguess!(t̄, q̄, t, q, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    sol = (t = t, q = q)
    hist = (t = [t̄], q = [q̄])
    extrapolate!(sol, hist, problem, extrap)
    (t = t, q = q)
end

function initialguess!(t, q, solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    initialguess!(previous(solstep)..., t, q, problem, extrap)
end



function initialguess!(t̄, q̄, p̄, t, q, p, q̇, ṗ, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    sol = (t = t, q = q, p = p)
    hist = (t = [t̄], q = [q̄], p = [p̄])
    extrapolate!(sol, hist, problem, extrap)
    functions(problem).v(q̇, t, q, p, parameters(problem))
    functions(problem).f(ṗ, t, q, p, parameters(problem))
    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t̄, q̄, p̄, t, q, p, q̇, ṗ, problem::AbstractProblemIODE, extrap::MidpointExtrapolation)
    sol = (t = t, q = q, p = p)
    hist = (t = [t̄], q = [q̄], p = [p̄])
    extrapolate!(sol, hist, problem, extrap)
    functions(problem).v̄(q̇, t, q, p, parameters(problem))
    functions(problem).f̄(ṗ, t, q, q̇, parameters(problem))
    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepPODE, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    initialguess!(previous(solstep)..., t, q, p, q̇, ṗ, problem, extrap)
end

function initialguess!(t̄, q̄, p̄, t, q, p, problem::Union{AbstractProblemPODE,AbstractProblemIODE}, extrap::MidpointExtrapolation)
    sol = (t = t, q = q, p = p)
    hist = (t = [t̄], q = [q̄], p = [p̄])
    extrapolate!(sol, hist, problem, extrap)
    (t = t, q = q, p = p)
end

function initialguess!(t, q, p, solstep::SolutionStepPODE, problem::Union{AbstractProblemPODE,AbstractProblemIODE}, extrap::MidpointExtrapolation)
    initialguess!(previous(solstep)..., t, q, p, problem, extrap)
end

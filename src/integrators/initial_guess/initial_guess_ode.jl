
function initialguess!(t, q, q̇, solstep::SolutionStepODE, ::AbstractProblemODE, ::NoInitialGuess)
    t  = solstep.t̄[1]
    q .= solstep.q̄[1]
    q̇ .= solstep.v̄[1]
    (t = t, q = q, v = q̇)
end


function initialguess!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, iguess::HermiteExtrapolation)
    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
        q̇ .= q̇₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, iguess)
    end
    (t = t, q = q, v = q̇)
end

function initialguess!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, iguess::HermiteExtrapolation)
    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, iguess)
    end
    (t = t, q = q)
end

function initialguess!(t, q, q̇, solstep::SolutionStepODE, ::AbstractProblemODE, extrap::HermiteExtrapolation)
    initialguess!(solstep.t̄[2], solstep.q̄[2], solstep.v̄[2], solstep.t̄[1], solstep.q̄[1], solstep.v̄[1], t, q, q̇, extrap)
end

function initialguess!(t, q, solstep::SolutionStepODE, ::AbstractProblemODE, extrap::HermiteExtrapolation)
    initialguess!(solstep.t̄[2], solstep.q̄[2], solstep.v̄[2], solstep.t̄[1], solstep.q̄[1], solstep.v̄[1], t, q, extrap)
end


function initialguess!(t̄, q̄, t, q, q̇, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, t, q, problem, extrap)
    functions(problem).v(q̇, t, q)
    (t = t, q = q, v = q̇)
end

function initialguess!(t, q, q̇, solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    initialguess!(solstep.t̄[1], solstep.q̄[1], t, q, q̇, problem, extrap)
end

function initialguess!(t̄, q̄, t, q, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, t, q, problem, extrap)
    (t = t, q = q)
end

function initialguess!(t, q, solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::MidpointExtrapolation)
    initialguess!(solstep.t̄[1], solstep.q̄[1], t, q, problem, extrap)
end

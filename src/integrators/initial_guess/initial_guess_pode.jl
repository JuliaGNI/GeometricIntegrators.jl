
function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepPODE, ::AbstractProblemPODE, ::NoInitialGuess)
    t  = solstep.t̄[1]
    q .= solstep.q̄[1]
    p .= solstep.p̄[1]
    q̇ .= solstep.v̄[1]
    ṗ .= solstep.f̄[1]
    (t, q, p, q̇, ṗ)
end


function initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, q̇, ṗ, iguess::HermiteExtrapolation)
    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
        q̇ .= q̇₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, iguess)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical!"
        p .= p₁
        ṗ .= ṗ₁
    else
        extrapolate!(t₀, p₀, ṗ₀, t₁, p₁, ṗ₁, t, p, ṗ, iguess)
    end

    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, iguess::HermiteExtrapolation)
    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, iguess)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical!"
        p .= p₁
    else
        extrapolate!(t₀, p₀, ṗ₀, t₁, p₁, ṗ₁, t, p, iguess)
    end

    (t = t, q = q, p = p)
end

function initialguess!(t, q, q̇, solstep::SolutionStepPODE, ::AbstractProblemPODE, extrap::HermiteExtrapolation)
    initialguess!(solstep.t̄[2], solstep.q̄[2], solstep.p̄[2], solstep.v̄[2], solstep.f̄[2], solstep.t̄[1], solstep.q̄[1], solstep.p̄[1], solstep.v̄[1], solstep.f̄[1], t, q, p, q̇, ṗ, extrap)
end

function initialguess!(t, q, solstep::SolutionStepPODE, ::AbstractProblemPODE, extrap::HermiteExtrapolation)
    initialguess!(solstep.t̄[2], solstep.q̄[2], solstep.p̄[2], solstep.v̄[2], solstep.f̄[2], solstep.t̄[1], solstep.q̄[1], solstep.p̄[1], solstep.v̄[1], solstep.f̄[1], t, q, p, extrap)
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

function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepODE, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    initialguess!(solstep.t̄[1], solstep.q̄[1], solstep.p̄[1], t, q, p, q̇, ṗ, problem, extrap)
end

function initialguess!(t̄, q̄, p̄, t, q, p, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    extrapolate!(t̄, q̄, t, q, problem, extrap)
    (t = t, q = q, p = p)
end

function initialguess!(t, q, p, solstep::SolutionStepODE, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    initialguess!(solstep.t̄[1], solstep.q̄[1], solstep.p̄[1], t, q, p, problem, extrap)
end

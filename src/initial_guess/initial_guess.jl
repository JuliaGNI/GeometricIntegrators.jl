
abstract type InitialGuess end


struct NoInitialGuess <: InitialGuess end

function initialguess!(t, q, q̇, solstep::SolutionStepODE, ::AbstractProblemODE, ::NoInitialGuess)
    t  = solstep.t̄[1]
    q .= solstep.q̄[1]
    q̇ .= solstep.v̄[1]
    (t = t, q = q, v = q̇)
end

function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepPODE, ::AbstractProblemPODE, ::NoInitialGuess)
    t  = solstep.t̄[1]
    q .= solstep.q̄[1]
    p .= solstep.p̄[1]
    q̇ .= solstep.v̄[1]
    ṗ .= solstep.f̄[1]
    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t₀, q₀, p₀, λ₀, q̇₀, ṗ₀, u₀, g₀, t₁, q₁, p₁, λ₁, q̇₁, ṗ₁, u₁, g₁, t, q, p, iguess::Union{Extrapolation,InitialGuess}; kwargs...)
    initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, iguess; kwargs...)
end

function initialguess!(t₀, q₀, p₀, λ₀, q̇₀, ṗ₀, u₀, g₀, t₁, q₁, p₁, λ₁, q̇₁, ṗ₁, u₁, g₁, t, q, p, q̇, ṗ, iguess::Union{Extrapolation,InitialGuess}; kwargs...)
    initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, q̇, ṗ, iguess; kwargs...)
end

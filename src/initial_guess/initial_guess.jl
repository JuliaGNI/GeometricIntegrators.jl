
abstract type InitialGuess end


struct NoInitialGuess <: InitialGuess end

function initialguess!(t, q, q̇, solstep::SolutionStepODE, ::AbstractProblemODE, ::NoInitialGuess)
    t  = solstep.t̄
    q .= solstep.q̄
    q̇ .= solstep.v̄
    (t = t, q = q, v = q̇)
end

function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepPODE, ::Union{AbstractProblemPODE,AbstractProblemIODE}, ::NoInitialGuess)
    t  = solstep.t̄
    q .= solstep.q̄
    p .= solstep.p̄
    q̇ .= solstep.v̄
    ṗ .= solstep.f̄
    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t₀, q₀, p₀, λ₀, μ₀, q̇₀, ṗ₀, u₀, g₀, t₁, q₁, p₁, λ₁, μ₁, q̇₁, ṗ₁, u₁, g₁, t, q, p, iguess::Union{Extrapolation,InitialGuess}; kwargs...)
    initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, iguess; kwargs...)
end

function initialguess!(t₀, q₀, p₀, λ₀, μ₀, q̇₀, ṗ₀, u₀, g₀, t₁, q₁, p₁, λ₁, μ₁, q̇₁, ṗ₁, u₁, g₁, t, q, p, q̇, ṗ, iguess::Union{Extrapolation,InitialGuess}; kwargs...)
    initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, q̇, ṗ, iguess; kwargs...)
end

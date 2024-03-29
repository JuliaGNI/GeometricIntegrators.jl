
function initialguess!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, iguess::HermiteExtrapolation; nowarn = false)
    if q₀ == q₁
        nowarn || @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
        q̇ .= q̇₁
    elseif t == t₁
        q .= q₁
        q̇ .= q̇₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, iguess)
    end
    (t = t, q = q, v = q̇)
end

function initialguess!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, iguess::HermiteExtrapolation; nowarn = false)
    if q₀ == q₁
        nowarn || @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
    elseif t == t₁
        q .= q₁
        q̇ .= q̇₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, iguess)
    end
    (t = t, q = q)
end

function initialguess!(t, q, q̇, solstep::Union{SolutionStepODE,SolutionStepDAE}, ::AbstractProblemODE, extrap::HermiteExtrapolation; kwargs...)
    t₀, q₀, q̇₀ = history(solstep, 2).t, history(solstep, 2).q, history(solstep, 2).v
    t₁, q₁, q̇₁ = history(solstep, 1).t, history(solstep, 1).q, history(solstep, 1).v
    initialguess!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, extrap; kwargs...)
end

function initialguess!(t, q, solstep::Union{SolutionStepODE,SolutionStepDAE}, ::AbstractProblemODE, extrap::HermiteExtrapolation; kwargs...)
    t₀, q₀, q̇₀ = history(solstep, 2).t, history(solstep, 2).q, history(solstep, 2).v
    t₁, q₁, q̇₁ = history(solstep, 1).t, history(solstep, 1).q, history(solstep, 1).v
    initialguess!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, extrap; kwargs...)
end



function initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, q̇, ṗ, iguess::HermiteExtrapolation; nowarn = false)
    if q₀ == q₁
        nowarn || @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
        q̇ .= q̇₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, q̇, iguess)
    end

    if p₀ == p₁
        nowarn || @warn "p₀ and p₁ in initial guess are identical!"
        p .= p₁
        ṗ .= ṗ₁
    else
        extrapolate!(t₀, p₀, ṗ₀, t₁, p₁, ṗ₁, t, p, ṗ, iguess)
    end

    (t = t, q = q, p = p, v = q̇, f = ṗ)
end

function initialguess!(t₀, q₀, p₀, q̇₀, ṗ₀, t₁, q₁, p₁, q̇₁, ṗ₁, t, q, p, iguess::HermiteExtrapolation; nowarn = false)
    if q₀ == q₁
        nowarn || @warn "q₀ and q₁ in initial guess are identical!"
        q .= q₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t, q, iguess)
    end

    if p₀ == p₁
        nowarn || @warn "p₀ and p₁ in initial guess are identical!"
        p .= p₁
    else
        extrapolate!(t₀, p₀, ṗ₀, t₁, p₁, ṗ₁, t, p, iguess)
    end

    (t = t, q = q, p = p)
end

function initialguess!(t, q, p, q̇, ṗ, solstep::Union{SolutionStepPODE,SolutionStepPDAE}, ::Union{AbstractProblemPODE,AbstractProblemIODE}, extrap::HermiteExtrapolation; kwargs...)
    initialguess!(history(solstep, 2)..., history(solstep, 1)..., t, q, p, q̇, ṗ, extrap; kwargs...)
end

function initialguess!(t, q, p, solstep::Union{SolutionStepPODE,SolutionStepPDAE}, ::Union{AbstractProblemPODE,AbstractProblemIODE}, extrap::HermiteExtrapolation; kwargs...)
    initialguess!(history(solstep, 2)..., history(solstep, 1)..., t, q, p, extrap; kwargs...)
end

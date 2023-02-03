
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

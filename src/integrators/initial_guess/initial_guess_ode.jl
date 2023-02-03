"""
`InitialGuessODE`: Initial guess for ordinary differential equations

### Fields

* `int`: interpolation structure
* `v`:   vector  field
* `Δt`:  time step
* `s`:   number of extrapolation stages (for initialisation)
"""
struct InitialGuessODE{TT, VT, IT <: Extrapolation}
    int::IT
    v::VT
    Δt::TT
    s::Int

    function InitialGuessODE(int::IT, v::VT, Δt::TT) where {TT <: Real, VT <: Function, IT <: Extrapolation}
        new{TT,VT,IT}(int, v, Δt, get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessODE(interp::Type{<:Extrapolation}, v::Function, Δt)
    int = interp()
    InitialGuessODE(int, v, Δt)
end

function InitialGuessODE(interp::Type{<:Extrapolation}, problem::Union{ODEProblem,DAEProblem,IODEProblem,LODEProblem}, Δt = tstep(problem))
    InitialGuessODE(interp, _get_v̄(equation(problem), parameters(problem)), Δt)
end

InitialGuess(interp, problem::Union{ODEProblem,DAEProblem}, Δt = tstep(problem)) = InitialGuessODE(interp, problem, Δt)
InitialGuess(problem::Union{ODEProblem,DAEProblem}, Δt = tstep(problem)) = InitialGuessODE(get_config(:ig_extrapolation), problem, Δt)

InitialGuess(interp, equation::ODEProblem, Δt = tstep(equation)) = InitialGuessODE(interp, equation, Δt)
InitialGuess(interp, equation::DAEProblem, Δt = tstep(equation)) = InitialGuessODE(interp, equation, Δt)


Base.:(==)(ig1::InitialGuessODE{TT1}, ig2::InitialGuessODE{TT2}) where {TT1,TT2}= (
                                TT1 == TT2
                             && ig1.int == ig2.int
                             && ig1.v   == ig2.v
                             && ig1.Δt  == ig2.Δt
                             && ig1.s   == ig2.s)


"Initialise initial guess, i.e., given t₀, t₁, q₁ compute q₀, v₀, v₁."
function initialize!(ig::InitialGuessODE{TT},
                t₁::TT,
                q₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                t₀::TT,
                q₀::SolutionVector{DT},
                v₀::SolutionVector{DT}) where {DT,TT}

    extrapolate_ode!(t₀, q₀, t₁, q₁, ig.v, MidpointExtrapolation(ig.s))

    ig.v(v₀, t₀, q₀)
    ig.v(v₁, t₁, q₁)
end


"compute vector field of new solution"
function Solutions.update_vector_fields!(ig::InitialGuessODE{TT}, t₁::TT,
                               q₁::SolutionVector{DT},
                               v₁::Vector{DT}) where {DT,TT}
    ig.v(v₁, t₁, q₁)
end

function GeometricBase.evaluate!(ig::InitialGuessODE{TT},
                q₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                q₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                guess::SolutionVector{DT},
                c::TT=one(TT)) where {DT,TT}

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁."
        guess .= q₁
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, (one(TT)+c)*ig.Δt, guess)
    end
end

function GeometricBase.evaluate!(ig::InitialGuessODE{TT},
                q₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                q₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                guess_q::SolutionVector{DT},
                guess_v::SolutionVector{DT},
                c::TT=one(TT)) where {DT,TT}
                @assert length(guess_q) == length(guess_v)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0."
        guess_q .= q₁
        guess_v .= 0
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, (one(TT)+c)*ig.Δt, guess_q, guess_v)
    end
end
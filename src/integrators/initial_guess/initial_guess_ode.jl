"""
`InitialGuessODE`: Initial guess for ordinary differential equations

### Fields

* `int`: interpolation structure
* `v`:   vector  field
* `Δt`:  time step
* `s`:   number of extrapolation stages (for initialisation)
"""
mutable struct InitialGuessODE{DT, TT, VT, IT <: Interpolator}
    int::IT
    v::VT
    Δt::TT
    s::Int

    function InitialGuessODE{DT,TT,VT,IT}(interp, v, Δt) where {DT,TT,VT,IT}
        new(interp, v, Δt, get_config(:ig_extrapolation_stages))
    end
end


function InitialGuessODE(interp, equation::ODE{DT,TT,VT}, Δt::TT) where {DT,TT,VT}
    int = interp(zero(DT), one(DT), Δt, ndims(equation))
    InitialGuessODE{DT, TT, VT, interp}(int, equation.v, Δt)
end

function InitialGuessODE(interp, equation::IODE{DT,TT,ΘT,FT,GT,VT}, Δt::TT) where {DT,TT,ΘT,FT,GT,VT}
    int = interp(zero(DT), one(DT), Δt, ndims(equation))
    InitialGuessODE{DT, TT, VT, interp}(int, equation.v, Δt)
end

function InitialGuessODE(interp, equation::VODE{DT,TT,AT,FT,GT,VT}, Δt::TT) where {DT,TT,AT,FT,GT,VT}
    int = interp(zero(DT), one(DT), Δt, ndims(equation))
    InitialGuessODE{DT, TT, VT, interp}(int, equation.v, Δt)
end


"Initialise initial guess, i.e., given t₀, t₁, q₁ compute q₀, v₀, v₁."
function initialize!(ig::InitialGuessODE{DT,TT,VT,IT},
                t₁::TT,
                q₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                t₀::TT,
                q₀::SolutionVector{DT},
                v₀::SolutionVector{DT}) where {DT,TT,VT,IT}
    midpoint_extrapolation(ig.v, t₁, t₀, q₁, q₀, ig.s)
    ig.v(t₀, q₀, v₀)
    ig.v(t₁, q₁, v₁)
end


"compute vector field of new solution"
function update!(ig::InitialGuessODE{DT,TT}, t₁::TT, q₁::SolutionVector{DT}, v₁::Vector{DT}) where {DT,TT}
    ig.v(t₁, q₁, v₁)
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT},
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
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c, guess)
    end
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT},
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
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c, guess_q, guess_v)
    end
end

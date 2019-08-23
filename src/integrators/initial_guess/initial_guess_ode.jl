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

function InitialGuessODE(interp, equation::VODE{DT,TT,AT,FT,GT,VT}, Δt::TT) where {DT,TT,AT,FT,GT,VT}
    int = interp(zero(DT), one(DT), Δt, ndims(equation))
    InitialGuessODE{DT, TT, VT, interp}(int, equation.v, Δt)
end


"Initialise initial guess, i.e., given t₀, t₁, q₁ compute q₀, v₀, v₁."
function initialize!(ig::InitialGuessODE{DT,TT,VT,IT},
                t₁::TT,
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                t₀::TT,
                q₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}}) where {DT,TT,VT,IT}
    midpoint_extrapolation(ig.v, t₁, t₀, q₁, q₀, ig.s)
    ig.v(t₀, q₀, v₀)
    ig.v(t₁, q₁, v₁)
end


"compute vector field of new solution"
function update!(ig::InitialGuessODE{DT,TT}, t₁::TT, q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, v₁::Vector{DT}) where {DT,TT}
    ig.v(t₁, q₁, v₁)
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT},
                q₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                c::TT=one(TT)) where {DT,TT}

    if q₀ == q₁
        tt₁ = zero(ig.Δt)
        tt₀ = - ig.Δt
        tq₀ = zero(q₀)
        tv₀ = zero(p₀)

        midpoint_extrapolation(ig.v, ig.f, tt₁, tt₀, q₁, tq₀, ig.s)
        ig.v(tt₀, tq₀, tp₀, tv₀)
        evaluate!(ig.int, tq₀, q₁, tv₀, v₁, one(TT)+c, guess)
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c, guess)
    end
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT},
                q₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_q::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_v::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                c::TT=one(TT)) where {DT,TT}
                @assert length(guess_q) == length(guess_v)

    if q₀ == q₁
        tt₁ = zero(ig.Δt)
        tt₀ = - ig.Δt # TODO # This is not the right time!
        tq₀ = zero(q₀)
        tv₀ = zero(p₀)

        midpoint_extrapolation(ig.v, ig.f, tt₁, tt₀, q₁, tq₀, ig.s)
        ig.v(tt₀, tq₀, tp₀, tv₀)
        evaluate!(ig.int, tq₀, q₁, tv₀, v₁, one(TT)+c, guess_q, guess_v)
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c, guess_q, guess_v)
    end
end

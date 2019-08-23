"""
`InitialGuessPODE`: Initial guess for partitioned ordinary differential equations

### Fields

* `int`: interpolation structure
* `v`:   vector field for q
* `f`:   vector field for p
* `Δt`:  time step
* `s`:   number of extrapolation stages (for initialisation)
"""
mutable struct InitialGuessPODE{DT, TT, VT, FT, IT <: Interpolator}
    int::IT
    v::VT
    f::FT
    Δt::TT
    s::Int

    function InitialGuessPODE{DT,TT,VT,FT,IT}(interp, v, f, Δt) where {DT,TT,VT,FT,IT}
        new(interp, v, f, Δt, get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessPODE(interp, equ::PODE{DT,TT,VT,FT}, Δt::TT) where {DT,TT,VT,FT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v, equ.f, Δt)
end

function InitialGuessPODE(interp, equ::IODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v, equ.f, Δt)
end

function InitialGuessPODE(interp, equ::VODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v, equ.f, Δt)
end

function InitialGuessPODE(interp, equ::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT}, Δt::TT) where {DT,TT,FT,PT,UT,GT,ϕT,VT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v, equ.f, Δt)
end


"Initialise initial guess, i.e., given t₀, t₁, q₁ compute q₀, v₀, v₁."
function initialize!(ig::InitialGuessPODE{DT,TT,VT,IT},
                t₁::TT,
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                f₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                t₀::TT,
                q₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                f₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}}) where {DT,TT,VT,IT}

    midpoint_extrapolation(ig.v, ig.f, t₁, t₀, q₁, q₀, p₁, p₀, ig.s)

    ig.v(t₀, q₀, p₀, v₀)
    ig.v(t₁, q₁, p₁, v₁)
    ig.f(t₀, q₀, v₀, f₀)
    ig.f(t₁, q₁, v₁, f₁)
end


function update!(ig::InitialGuessPODE{DT,TT}, t₁::TT,
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Vector{DT},
                f₁::Vector{DT}) where {DT,TT}
    ig.v(t₁, q₁, p₁, v₁)
    ig.f(t₁, q₁, v₁, f₁)
end


function CommonFunctions.evaluate!(ig::InitialGuessPODE{DT,TT},
                q₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                f₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                f₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_q::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_p::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                c_q::TT=one(TT), c_p::TT=one(TT)) where {DT,TT}
                @assert length(guess_q) == length(guess_p)

    if q₀ == q₁
        tt₁ = zero(ig.Δt)
        tt₀ = - ig.Δt # TODO # This is not the right time!
        tq₀ = zero(q₀)
        tp₀ = zero(p₀)
        tv₀ = zero(p₀)

        midpoint_extrapolation(ig.v, ig.f, tt₁, tt₀, q₁, tq₀, p₁, tp₀, ig.s)
        ig.v(tt₀, tq₀, tp₀, tv₀)
        evaluate!(ig.int, tq₀, q₁, tv₀, v₁, one(TT)+c_q, guess_q)
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q)
    end

    if p₀ == p₁
        tt₁ = zero(ig.Δt)
        tt₀ = - ig.Δt # TODO # This is not the right time!
        tq₀ = zero(q₀)
        tp₀ = zero(p₀)
        tv₀ = zero(p₀)
        tf₀ = zero(p₀)

        midpoint_extrapolation(ig.v, ig.f, tt₁, tt₀, q₁, tq₀, p₁, tp₀, ig.s)
        ig.v(tt₀, tq₀, tp₀, tv₀)
        ig.f(tt₀, tq₀, tv₀, tf₀)
        evaluate!(ig.int, tp₀, p₁, tf₀, f₁, one(TT)+c_p, guess_p)
    else
        evaluate!(ig.int, p₀, p₁, f₀, f₁, one(TT)+c_p, guess_p)
    end
end


function CommonFunctions.evaluate!(ig::InitialGuessPODE{DT,TT},
                q₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                f₀::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                p₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                v₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                f₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_q::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_p::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_v::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                guess_f::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                c_q::TT=one(TT), c_p::TT=one(TT)) where {DT,TT}
                @assert length(guess_q) == length(guess_v)
                @assert length(guess_p) == length(guess_f)
                @assert length(guess_q) == length(guess_p)

    if q₀ == q₁
        tt₁ = zero(ig.Δt)
        tt₀ = - ig.Δt # TODO # This is not the right time!
        tq₀ = zero(q₀)
        tp₀ = zero(p₀)
        tv₀ = zero(p₀)

        midpoint_extrapolation(ig.v, ig.f, tt₁, tt₀, q₁, tq₀, p₁, tp₀, ig.s)
        ig.v(tt₀, tq₀, tp₀, tv₀)
        evaluate!(ig.int, tq₀, q₁, tv₀, v₁, one(TT)+c_q, guess_q, guess_v)
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q, guess_v)
    end

    if p₀ == p₁
        tt₁ = zero(ig.Δt)
        tt₀ = - ig.Δt # TODO # This is not the right time!
        tq₀ = zero(q₀)
        tp₀ = zero(p₀)
        tv₀ = zero(p₀)
        tf₀ = zero(p₀)

        midpoint_extrapolation(ig.v, ig.f, tt₁, tt₀, q₁, tq₀, p₁, tp₀, ig.s)
        ig.v(tt₀, tq₀, tp₀, tv₀)
        ig.f(tt₀, tq₀, tv₀, tf₀)
        evaluate!(ig.int, tp₀, p₁, tf₀, f₁, one(TT)+c_p, guess_p, guess_f)
    else
        evaluate!(ig.int, p₀, p₁, f₀, f₁, one(TT)+c_p, guess_p, guess_f)
    end
end

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
                c_q::TT) where {DT,TT}

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁."
        guess .= q₁
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q)
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
                guess_v::Union{Vector{DT}, Vector{TwicePrecision{DT}}},
                c_q::TT) where {DT,TT}
    @assert length(guess_q) == length(guess_v)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0."
        guess_q .= q₁
        guess_v .= 0
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q, guess_v)
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
                c_q::TT, c_p::TT) where {DT,TT}
    @assert length(guess_q) == length(guess_p)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁."
        guess_q .= q₁
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical! Setting p=p₁."
        guess_p .= p₁
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
                c_q::TT, c_p::TT) where {DT,TT}
    @assert length(guess_q) == length(guess_v)
    @assert length(guess_p) == length(guess_f)
    @assert length(guess_q) == length(guess_p)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0."
        guess_q .= q₁
        guess_v .= 0
    else
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q, guess_v)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical! Setting p=p₁ and f=0."
        guess_p .= p₁
        guess_f .= 0
    else
        evaluate!(ig.int, p₀, p₁, f₀, f₁, one(TT)+c_p, guess_p, guess_f)
    end
end

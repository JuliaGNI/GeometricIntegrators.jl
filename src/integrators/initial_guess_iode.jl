
using ..CommonFunctions
using ..Interpolation
using ..Tableaus


type InitialGuessIODE{DT, TT, VT, FT, IT <: Interpolator}
    int::IT

    v::VT
    f::FT

    Δt::TT
    t₀::TT
    t₁::TT

    q₀::Vector{DT}
    q₁::Vector{DT}
    v₀::Vector{DT}
    v₁::Vector{DT}
    p₀::Vector{DT}
    p₁::Vector{DT}
    f₀::Vector{DT}
    f₁::Vector{DT}
end

function InitialGuessIODE{DT,TT,ΑT,FT,GT,VT}(interp, equ::IODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT)
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, 0, 0,
                                         zeros(DT, equ.d), zeros(DT, equ.d),
                                         zeros(DT, equ.d), zeros(DT, equ.d),
                                         zeros(DT, equ.d), zeros(DT, equ.d),
                                         zeros(DT, equ.d), zeros(DT, equ.d))
end

function InitialGuessIODE{DT,TT,FT,PT,UT,GT,ϕT,VT}(interp, equ::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT}, Δt::TT)
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, 0, 0,
                                         zeros(DT, equ.d), zeros(DT, equ.d),
                                         zeros(DT, equ.d), zeros(DT, equ.d),
                                         zeros(DT, equ.d), zeros(DT, equ.d),
                                         zeros(DT, equ.d), zeros(DT, equ.d))
end

function initialize!{DT,TT,VT,FT,IT}(ig::InitialGuessIODE{DT,TT,VT,FT,IT},
                                     t₁::TT, q₁::Vector{DT}, p₁::Vector{DT})
    ig.t₁ = t₁
    # TODO Replace copy with EPRK4 step.
    simd_copy!(q₁, ig.q₁)
    simd_copy!(p₁, ig.p₁)

    ig.v(t₁-ig.Δt, ig.q₁, ig.p₁, ig.v₁)
    ig.f(t₁-ig.Δt, ig.q₁, ig.v₁, ig.f₁)
end

function update!{DT,TT,VT,FT,IT}(ig::InitialGuessIODE{DT,TT,VT,FT,IT},
                                 t₁::TT, q₁::Vector{DT}, p₁::Vector{DT})
    ig.t₀ = ig.t₁
    ig.t₁ = t₁

    simd_copy!(ig.q₁, ig.q₀)
    simd_copy!(ig.v₁, ig.v₀)
    simd_copy!(ig.p₁, ig.p₀)
    simd_copy!(ig.f₁, ig.f₀)

    simd_copy!(q₁, ig.q₁)
    simd_copy!(p₁, ig.p₁)

    ig.v(ig.t₁, ig.q₁, ig.p₁, ig.v₁)
    ig.f(ig.t₁, ig.q₁, ig.v₁, ig.f₁)
end

function CommonFunctions.evaluate!{DT,TT,VT,FT,IT}(ig::InitialGuessIODE{DT,TT,VT,FT,IT},
           guess_q::Vector{DT}, guess_p::Vector{DT}, guess_v::Vector{DT},
           c_q::TT=one(TT), c_p::TT=one(TT))

    @assert length(guess_q) == length(guess_p) == length(guess_v)

    evaluate!(ig.int, ig.q₀, ig.q₁, ig.v₀, ig.v₁, one(TT)+c_q, guess_q, guess_v)

    if ig.q₀ == ig.q₁
        guess_v .= 0
    end

    evaluate!(ig.int, ig.p₀, ig.p₁, ig.f₀, ig.f₁, one(TT)+c_p, guess_p)

    nothing
end

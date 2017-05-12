
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

    periodicity::Vector{DT}

    q₀::Vector{DT}
    q₁::Vector{DT}
    v₀::Vector{DT}
    v₁::Vector{DT}
    p₀::Vector{DT}
    p₁::Vector{DT}
    f₀::Vector{DT}
    f₁::Vector{DT}

    function InitialGuessIODE(interp, v, f, Δt, d, periodicity)
        if !(length(periodicity) == d)
            periodicity = zeros(DT, d)
        end
        new(interp, v, f, Δt, 0, 0, periodicity,
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d))
    end
end

function InitialGuessIODE{DT,TT,ΑT,FT,GT,VT}(interp, equ::IODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT; periodicity=[])
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.d, periodicity)
end

function InitialGuessIODE{DT,TT,FT,PT,UT,GT,ϕT,VT}(interp, equ::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT}, Δt::TT; periodicity=[])
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.d, periodicity)
end

function initialize!{DT,TT,VT,FT,IT}(ig::InitialGuessIODE{DT,TT,VT,FT,IT},
                                     t₁::TT, q₁::Vector{DT}, p₁::Vector{DT})
    ig.t₁ = t₁
    # TODO Replace copy with EPRK4 step.
    simd_copy!(q₁, ig.q₁)
    simd_copy!(p₁, ig.p₁)

    ig.v(t₁-ig.Δt, ig.q₁, ig.p₁, ig.v₁)
    ig.f(t₁-ig.Δt, ig.q₁, ig.v₁, ig.f₁)

    update!(ig, t₁, q₁, p₁)
end

function update!{DT,TT,VT,FT,IT}(ig::InitialGuessIODE{DT,TT,VT,FT,IT},
                                 t₁::TT, q₁::Vector{DT}, p₁::Vector{DT})
    local Δq::DT

    ig.t₀ = ig.t₁
    ig.t₁ = t₁

    simd_copy!(ig.q₁, ig.q₀)
    simd_copy!(ig.v₁, ig.v₀)
    simd_copy!(ig.p₁, ig.p₀)
    simd_copy!(ig.f₁, ig.f₀)

    simd_copy!(q₁, ig.q₁)
    simd_copy!(p₁, ig.p₁)

    # take care of periodic solutions
    for k in 1:length(ig.q₁)
        if ig.periodicity[k] ≠ 0
            Δq = ig.q₁[k] - ig.q₀[k]
            if ig.q₁[k] < 0 || ig.q₁[k] ≥ ig.periodicity[k]
                ig.q₁[k] -= fld(ig.q₁[k], ig.periodicity[k]) * ig.periodicity[k]
            end
            ig.q₀[k] = ig.q₁[k] - Δq
        end
    end

    # compute vector fields
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

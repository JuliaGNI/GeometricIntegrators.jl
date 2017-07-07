
using ..CommonFunctions
using ..Interpolation
using ..Tableaus


"""
Midpoint extrapolation method with arbitrary order p.

    v:  function to compute vector field
    f:  function to compute force  field
    t₀: initial time
    t₁: final   time
    q₀: initial positions
    p₀: initial momenta
    q₁: final   positions
    p₁: final   momenta
    s:  number of interpolations (order p=2s+2)
"""
function midpoint_extrapolation(v::Function, f::Function, t₀::TT, t₁::TT, q₀::Vector{DT}, p₀::Vector{DT}, q₁::Vector{DT}, p₁::Vector{DT}, s::Int) where {DT,TT}
    @assert size(q₀) == size(q₁) == size(p₀) == size(p₁)

    local F   = 2*collect(1:(s+1))
    local Δt  = t₁ - t₀
    local σ   = Δt ./ F
    local σ2  = σ.^2

    local qts = zeros(length(q₀), s+1)
    local pts = zeros(length(p₀), s+1)

    local qᵢ₁= zeros(q₀)
    local qᵢ₂= zeros(q₀)
    local qᵢₜ= zeros(q₀)

    local pᵢ₁= zeros(p₀)
    local pᵢ₂= zeros(p₀)
    local pᵢₜ= zeros(p₀)

    local v₀ = zeros(q₀)
    local vᵢ = zeros(q₀)

    local f₀ = zeros(p₀)
    local fᵢ = zeros(p₀)

    v(t₀, q₀, p₀, v₀)
    f(t₀, q₀, v₀, f₀)

    for i in 1:(s+1)
        tᵢ   = t₀ + σ[i]
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ + σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ + σ[i] .* f₀
        for j in 1:(F[i]-1)
            v(tᵢ, qᵢ₂, pᵢ₂, vᵢ)
            f(tᵢ, qᵢ₂, vᵢ,  fᵢ)
            qᵢₜ .= qᵢ₁ + 2σ[i] .* vᵢ
            qᵢ₁ .= qᵢ₂
            qᵢ₂ .= qᵢₜ
            pᵢₜ .= pᵢ₁ + 2σ[i] .* fᵢ
            pᵢ₁ .= pᵢ₂
            pᵢ₂ .= pᵢₜ
        end
        for k in indices(qts,1)
            qts[k,i] += qᵢ₂[k]
        end
        for k in indices(pts,1)
            pts[k,i] += pᵢ₂[k]
        end
    end

    aitken_neville(σ2, qts, zero(TT), q₁)
    aitken_neville(σ2, pts, zero(TT), p₁)
end



mutable struct InitialGuessIODE{DT, TT, VT, FT, IT <: Interpolator}
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

    s::Int

    function InitialGuessIODE{DT,TT,VT,FT,IT}(interp, v, f, Δt, d, periodicity) where {DT,TT,VT,FT,IT}
        if !(length(periodicity) == d)
            periodicity = zeros(DT, d)
        end
        new(interp, v, f, Δt, 0, 0, periodicity,
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d),
            get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessIODE(interp, equ::IODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT; periodicity=[]) where {DT,TT,ΑT,FT,GT,VT}
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.d, periodicity)
end

function InitialGuessIODE(interp, equ::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT}, Δt::TT; periodicity=[]) where {DT,TT,FT,PT,UT,GT,ϕT,VT}
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.d, periodicity)
end


function initialize!(ig::InitialGuessIODE{DT,TT,VT,FT,IT}, t₁::TT, q₁::Vector{DT}, p₁::Vector{DT}) where {DT,TT,VT,FT,IT}
    ig.t₀  = t₁ - ig.Δt
    ig.t₁  = t₁
    ig.q₁ .= q₁
    ig.p₁ .= p₁

    midpoint_extrapolation(ig.v, ig.f, ig.t₁, ig.t₀, ig.q₁, ig.p₁, ig.q₀, ig.p₀, ig.s)

    ig.v(ig.t₀, ig.q₀, ig.p₀, ig.v₀)
    ig.v(ig.t₁, ig.q₁, ig.p₁, ig.v₁)
    ig.f(ig.t₀, ig.q₀, ig.v₀, ig.f₀)
    ig.f(ig.t₁, ig.q₁, ig.v₁, ig.f₁)
end

function update!(ig::InitialGuessIODE{DT,TT,VT,FT,IT}, t₁::TT, q₁::Vector{DT}, p₁::Vector{DT}) where {DT,TT,VT,FT,IT}
    local Δq::DT

    ig.t₀ = ig.t₁
    ig.t₁ = t₁

    ig.q₀ .= ig.q₁
    ig.v₀ .= ig.v₁
    ig.p₀ .= ig.p₁
    ig.f₀ .= ig.f₁

    ig.q₁ .= q₁
    ig.p₁ .= p₁

    # compute vector fields
    ig.v(ig.t₁, ig.q₁, ig.p₁, ig.v₁)
    ig.f(ig.t₁, ig.q₁, ig.v₁, ig.f₁)

    # take care of periodic solutions
    for k in eachindex(ig.q₁)
        if ig.periodicity[k] ≠ 0
            Δq = ig.q₁[k] - ig.q₀[k]
            # if ig.q₁[k] < 0 || ig.q₁[k] ≥ ig.periodicity[k]
            #     ig.q₁[k] -= fld(ig.q₁[k], ig.periodicity[k]) * ig.periodicity[k]
            # end
            while ig.q₁[k] < 0
                ig.q₁[k] += ig.periodicity[k]
            end
            while ig.q₁[k] ≥ ig.periodicity[k]
                ig.q₁[k] -= ig.periodicity[k]
            end
            ig.q₀[k] = ig.q₁[k] - Δq
        end
    end
end

function CommonFunctions.evaluate!(ig::InitialGuessIODE{DT,TT,VT,FT,IT},
           guess_q::Vector{DT}, guess_p::Vector{DT}, guess_v::Vector{DT},
           c_q::TT=one(TT), c_p::TT=one(TT)) where {DT,TT,VT,FT,IT}

    @assert length(guess_q) == length(guess_p) == length(guess_v)

    evaluate!(ig.int, ig.q₀, ig.q₁, ig.v₀, ig.v₁, one(TT)+c_q, guess_q, guess_v)

    if ig.q₀ == ig.q₁
        warn("q₀ and q₁ in initial guess are identital!")
        guess_v .= 0
    end

    evaluate!(ig.int, ig.p₀, ig.p₁, ig.f₀, ig.f₁, one(TT)+c_p, guess_p)

    nothing
end

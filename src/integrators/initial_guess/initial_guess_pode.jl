
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

    local qts = zeros(eltype(q₀), length(q₀), s+1)
    local pts = zeros(eltype(p₀), length(p₀), s+1)

    local qᵢ₁= zero(q₀)
    local qᵢ₂= zero(q₀)
    local qᵢₜ= zero(q₀)

    local pᵢ₁= zero(p₀)
    local pᵢ₂= zero(p₀)
    local pᵢₜ= zero(p₀)

    local v₀ = zero(q₀)
    local vᵢ = zero(q₀)

    local f₀ = zero(p₀)
    local fᵢ = zero(p₀)

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
        for k in axes(qts,1)
            qts[k,i] += qᵢ₂[k]
        end
        for k in axes(pts,1)
            pts[k,i] += pᵢ₂[k]
        end
    end

    aitken_neville(σ2, qts, zero(TT), q₁)
    aitken_neville(σ2, pts, zero(TT), p₁)
end



mutable struct InitialGuessPODE{DT, TT, VT, FT, IT <: Interpolator}
    int::IT

    v::VT
    f::FT

    Δt::TT

    periodicity::Vector{DT}

    t₀::Vector{TT}
    q₀::Vector{Vector{DT}}
    v₀::Vector{Vector{DT}}
    p₀::Vector{Vector{DT}}
    f₀::Vector{Vector{DT}}

    t₁::Vector{TT}
    q₁::Vector{Vector{DT}}
    v₁::Vector{Vector{DT}}
    p₁::Vector{Vector{DT}}
    f₁::Vector{Vector{DT}}

    s::Int

    function InitialGuessPODE{DT,TT,VT,FT,IT}(interp, v, f, Δt, m, d, periodicity) where {DT,TT,VT,FT,IT}
        if !(length(periodicity) == d)
            periodicity = zeros(DT, d)
        end

        t₀ = zeros(DT,m)
        q₀ = Array{Vector{DT}}(undef, m)
        v₀ = Array{Vector{DT}}(undef, m)
        p₀ = Array{Vector{DT}}(undef, m)
        f₀ = Array{Vector{DT}}(undef, m)

        t₁ = zeros(DT,m)
        q₁ = Array{Vector{DT}}(undef, m)
        v₁ = Array{Vector{DT}}(undef, m)
        p₁ = Array{Vector{DT}}(undef, m)
        f₁ = Array{Vector{DT}}(undef, m)

        for i in 1:m
            q₀[i] = zeros(DT,d)
            v₀[i] = zeros(DT,d)
            p₀[i] = zeros(DT,d)
            f₀[i] = zeros(DT,d)
            q₁[i] = zeros(DT,d)
            v₁[i] = zeros(DT,d)
            p₁[i] = zeros(DT,d)
            f₁[i] = zeros(DT,d)
        end


        new(interp, v, f, Δt, periodicity, t₀, q₀, v₀, p₀, f₀, t₁, q₁, v₁, p₁, f₁,
            get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessPODE(interp, equ::PODE{DT,TT,VT,FT}, Δt::TT) where {DT,TT,VT,FT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.n, equ.d, equ.periodicity)
end

function InitialGuessPODE(interp, equ::IODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.n, equ.d, equ.periodicity)
end

function InitialGuessPODE(interp, equ::VODE{DT,TT,ΑT,FT,GT,VT}, Δt::TT) where {DT,TT,ΑT,FT,GT,VT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.n, equ.d, equ.periodicity)
end

function InitialGuessPODE(interp, equ::IDAE{DT,TT,FT,PT,UT,GT,ϕT,VT}, Δt::TT) where {DT,TT,FT,PT,UT,GT,ϕT,VT}
    InitialGuessPODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, equ.d),
                                         equ.v, equ.f, Δt, equ.n, equ.d, equ.periodicity)
end


function initialize!(ig::InitialGuessPODE{DT,TT}, m::Int, t₁::TT, q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}}) where {DT,TT}
    ig.t₀[m]  = t₁ - ig.Δt
    ig.t₁[m]  = t₁
    ig.q₁[m] .= q₁
    ig.p₁[m] .= p₁

    midpoint_extrapolation(ig.v, ig.f, ig.t₁[m], ig.t₀[m], ig.q₁[m], ig.p₁[m], ig.q₀[m], ig.p₀[m], ig.s)

    ig.v(ig.t₀[m], ig.q₀[m], ig.p₀[m], ig.v₀[m])
    ig.v(ig.t₁[m], ig.q₁[m], ig.p₁[m], ig.v₁[m])
    ig.f(ig.t₀[m], ig.q₀[m], ig.v₀[m], ig.f₀[m])
    ig.f(ig.t₁[m], ig.q₁[m], ig.v₁[m], ig.f₁[m])
end


function update!(ig::InitialGuessPODE{DT,TT}, m::Int, t₁::TT, q₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p₁::Union{Vector{DT}, Vector{TwicePrecision{DT}}}) where {DT,TT}
    local Δq::DT

    ig.t₀[m] = ig.t₁[m]
    ig.t₁[m] = t₁

    ig.q₀[m] .= ig.q₁[m]
    ig.v₀[m] .= ig.v₁[m]
    ig.p₀[m] .= ig.p₁[m]
    ig.f₀[m] .= ig.f₁[m]

    ig.q₁[m] .= q₁
    ig.p₁[m] .= p₁

    # compute vector fields
    ig.v(ig.t₁[m], ig.q₁[m], ig.p₁[m], ig.v₁[m])
    ig.f(ig.t₁[m], ig.q₁[m], ig.v₁[m], ig.f₁[m])

    # take care of periodic solutions
    for k in eachindex(ig.q₁[m], ig.periodicity)
        if ig.periodicity[k] ≠ 0
            Δq = ig.q₁[m][k] - ig.q₀[m][k]
            # if ig.q₁[m][k] < 0 || ig.q₁[m][k] ≥ ig.periodicity[k]
            #     ig.q₁[m][k] -= fld(ig.q₁[m][k], ig.periodicity[k]) * ig.periodicity[k]
            # end
            while ig.q₁[m][k] < 0
                ig.q₁[m][k] += ig.periodicity[k]
            end
            while ig.q₁[m][k] ≥ ig.periodicity[k]
                ig.q₁[m][k] -= ig.periodicity[k]
            end
            ig.q₀[m][k] = ig.q₁[m][k] - Δq
        end
    end
end


function CommonFunctions.evaluate!(ig::InitialGuessPODE{DT,TT,VT,FT,IT}, m::Int,
           guess_q::Vector{DT}, guess_p::Vector{DT}, guess_v::Vector{DT},
           c_q::TT=one(TT), c_p::TT=one(TT)) where {DT,TT,VT,FT,IT}

    @assert length(guess_q) == length(guess_p) == length(guess_v)

    if ig.q₀[m] == ig.q₁[m]
        warn("q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0.")
        guess_q .= ig.q₁[m]
        guess_v .= 0
    else
        evaluate!(ig.int, ig.q₀[m], ig.q₁[m], ig.v₀[m], ig.v₁[m], one(TT)+c_q, guess_q, guess_v)
    end

    if ig.p₀[m] == ig.p₁[m]
        warn("p₀ and p₁ in initial guess are identical! Setting p=p₁.")
        guess_p .= ig.p₁[m]
    else
        evaluate!(ig.int, ig.p₀[m], ig.p₁[m], ig.f₀[m], ig.f₁[m], one(TT)+c_p, guess_p)
    end
end


function CommonFunctions.evaluate!(ig::InitialGuessPODE{DT,TT,VT,FT,IT}, m::Int,
           guess_q::Vector{DT}, guess_p::Vector{DT}, guess_v::Vector{DT}, guess_f::Vector{DT},
           c_q::TT=one(TT), c_p::TT=one(TT)) where {DT,TT,VT,FT,IT}

    @assert length(guess_q) == length(guess_p) == length(guess_v) == length(guess_f)

    if ig.q₀[m] == ig.q₁[m]
        warn("q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0.")
        guess_q .= ig.q₁[m]
        guess_v .= 0
    else
        evaluate!(ig.int, ig.q₀[m], ig.q₁[m], ig.v₀[m], ig.v₁[m], one(TT)+c_q, guess_q, guess_v)
    end

    if ig.p₀[m] == ig.p₁[m]
        warn("p₀ and p₁ in initial guess are identical! Setting p=p₁.")
        guess_p .= ig.p₁[m]
        guess_f .= 0
    else
        evaluate!(ig.int, ig.p₀[m], ig.p₁[m], ig.f₀[m], ig.f₁[m], one(TT)+c_p, guess_p, guess_f)
    end
end

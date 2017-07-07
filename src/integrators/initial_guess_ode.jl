
using ..CommonFunctions
using ..Interpolation


mutable struct InitialGuessODE{DT, TT, VT, IT <: Interpolator}
    int::IT

    v::VT

    Δt::TT
    t₀::TT
    t₁::TT

    periodicity::Vector{DT}

    q₀::Vector{DT}
    v₀::Vector{DT}
    q₁::Vector{DT}
    v₁::Vector{DT}

    s::Int

    function InitialGuessODE{DT,TT,VT,IT}(interp, v, Δt, d, periodicity) where {DT,TT,VT,IT}
        new(interp, v, Δt, 0, 0, periodicity,
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d),
            get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessODE(interp, equation::ODE{DT,TT,VT,N}, Δt::TT; periodicity=[]) where {DT,TT,VT,N}
    if N > 1
        equ = similar(equation, equation.q₀[:,1])
    else
        equ = equation
    end

    if !(length(periodicity) == equation.d)
        periodicity = zeros(DT, equation.d)
    end

    int = interp(zero(DT), one(DT), Δt, equ.d)
    InitialGuessODE{DT, TT, VT, interp}(int, equ.v, Δt, equ.d, periodicity)
end

function initialize!(ig::InitialGuessODE{DT,TT,VT,IT}, t₁::TT, q₁::Union{Vector{DT}, Vector{Double{DT}}}) where {DT,TT,VT,IT}
    ig.t₀ = t₁ - ig.Δt
    ig.t₁ = t₁
    ig.q₁ = q₁

    # euler_extrapolation(ig.v, ig.t₁, ig.t₀, ig.q₁, ig.q₀, ig.s)
    midpoint_extrapolation(ig.v, ig.t₁, ig.t₀, ig.q₁, ig.q₀, ig.s)

    ig.v(ig.t₀, ig.q₀, ig.v₀)
    ig.v(ig.t₁, ig.q₁, ig.v₁)
end

function update!(ig::InitialGuessODE{DT,TT,VT,IT}, t₁::TT, q₁::Union{Vector{DT}, Vector{Double{DT}}}) where {DT,TT,VT,IT}
    local Δq::DT

    ig.t₀ = ig.t₁
    ig.t₁ = t₁

    ig.q₀ .= ig.q₁
    ig.v₀ .= ig.v₁
    ig.q₁ .= q₁

    # take care of periodic solutions
    for k in eachindex(ig.q₁)
        if ig.periodicity[k] ≠ 0
            Δq = ig.q₁[k] - ig.q₀[k]
            if ig.q₁[k] < 0 || ig.q₁[k] ≥ ig.periodicity[k]
                ig.q₁[k] -= fld(ig.q₁[k], ig.periodicity[k]) * ig.periodicity[k]
            end
            ig.q₀[k] = ig.q₁[k] - Δq
        end
    end

    # compute vector field
    ig.v(t₁, ig.q₁, ig.v₁)
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT,VT,IT}, guess::Union{Vector{DT}, Vector{Double{DT}}}, c::TT=one(TT)) where {DT,TT,VT,IT}
    evaluate!(ig.int, ig.q₀, ig.q₁, ig.v₀, ig.v₁, one(TT)+c, guess)
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT,VT,IT},
           guess_q::Union{Vector{DT}, Vector{Double{DT}}}, guess_v::Vector{DT}, c::TT=one(TT)) where {DT,TT,VT,IT}

    @assert length(guess_q) == length(guess_v)

    evaluate!(ig.int, ig.q₀, ig.q₁, ig.v₀, ig.v₁, one(TT)+c, guess_q, guess_v)

    if ig.q₀ == ig.q₁
        warn("q₀ and q₁ in initial guess are identital!")
        guess_v .= 0
    end

    nothing
end


using ..CommonFunctions
using ..Interpolation


mutable struct InitialGuessODE{DT, TT, VT, IT <: Interpolator}
    int::IT

    v::VT

    Δt::TT

    periodicity::Vector{DT}

    t₀::Vector{TT}
    q₀::Vector{Vector{DT}}
    v₀::Vector{Vector{DT}}

    t₁::Vector{TT}
    q₁::Vector{Vector{DT}}
    v₁::Vector{Vector{DT}}

    s::Int

    function InitialGuessODE{DT,TT,VT,IT}(interp, v, Δt, m, d, periodicity) where {DT,TT,VT,IT}

        t₀ = zeros(DT,m)
        q₀ = Array{Vector{DT}}(m)
        v₀ = Array{Vector{DT}}(m)

        t₁ = zeros(DT,m)
        q₁ = Array{Vector{DT}}(m)
        v₁ = Array{Vector{DT}}(m)

        for i in 1:m
            q₀[i] = zeros(DT,d)
            v₀[i] = zeros(DT,d)
            q₁[i] = zeros(DT,d)
            v₁[i] = zeros(DT,d)
        end

        new(interp, v, Δt, periodicity, t₀, q₀, v₀, t₁, q₁, v₁,
            get_config(:ig_extrapolation_stages))
    end
end


function InitialGuessODE(interp, equation::ODE{DT,TT,VT}, Δt::TT) where {DT,TT,VT}
    int = interp(zero(DT), one(DT), Δt, equation.d)
    InitialGuessODE{DT, TT, VT, interp}(int, equation.v, Δt, equation.n, equation.d, equation.periodicity)
end

function InitialGuessODE(interp, equation::VODE{DT,TT,AT,FT,GT,VT}, Δt::TT) where {DT,TT,AT,FT,GT,VT}
    int = interp(zero(DT), one(DT), Δt, equation.d)
    InitialGuessODE{DT, TT, VT, interp}(int, equation.v, Δt, equation.n, equation.d, equation.periodicity)
end


function initialize!(ig::InitialGuessODE{DT,TT,VT,IT}, m::Int, t₁::TT, q₁::Union{Vector{DT}, Vector{Double{DT}}}) where {DT,TT,VT,IT}
    ig.t₀[m]  = t₁ - ig.Δt
    ig.t₁[m]  = t₁
    ig.q₁[m] .= q₁

    midpoint_extrapolation(ig.v, ig.t₁[m], ig.t₀[m], ig.q₁[m], ig.q₀[m], ig.s)

    ig.v(ig.t₀[m], ig.q₀[m], ig.v₀[m])
    ig.v(ig.t₁[m], ig.q₁[m], ig.v₁[m])
end


function update!(ig::InitialGuessODE{DT,TT,VT,IT}, m::Int, t₁::TT, q₁::Union{Vector{DT}, Vector{Double{DT}}}) where {DT,TT,VT,IT}
    local Δq::DT

    ig.t₀[m] = ig.t₁[m]
    ig.t₁[m] = t₁

    ig.q₀[m] .= ig.q₁[m]
    ig.v₀[m] .= ig.v₁[m]
    ig.q₁[m] .= q₁

    # compute vector field
    ig.v(ig.t₁[m], ig.q₁[m], ig.v₁[m])

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

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT,VT,IT}, m::Int, guess::Union{Vector{DT}, Vector{Double{DT}}}, c::TT=one(TT)) where {DT,TT,VT,IT}
    evaluate!(ig.int, ig.q₀[m], ig.q₁[m], ig.v₀[m], ig.v₁[m], one(TT)+c, guess)
end

function CommonFunctions.evaluate!(ig::InitialGuessODE{DT,TT,VT,IT}, m::Int,
           guess_q::Union{Vector{DT}, Vector{Double{DT}}}, guess_v::Vector{DT}, c::TT=one(TT)) where {DT,TT,VT,IT}

    @assert length(guess_q) == length(guess_v)

    if ig.q₀[m] == ig.q₁[m]
        warn("q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0.")
        guess_q .= ig.q₁[m]
        guess_v .= 0
    else
        evaluate!(ig.int, ig.q₀[m], ig.q₁[m], ig.v₀[m], ig.v₁[m], one(TT)+c, guess_q, guess_v)
    end
end

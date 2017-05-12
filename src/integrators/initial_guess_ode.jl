
using ..CommonFunctions
using ..Interpolation
using ..Tableaus


type InitialGuessODE{DT, TT, VT, IT <: Interpolator}
    int::IT

    rk4::IntegratorERK{DT,TT,VT}
    sol::SolutionODE{DT,TT,2}

    v::VT

    Δt::TT
    t₀::TT
    t₁::TT

    periodicity::Vector{DT}

    q₀::Vector{DT}
    v₀::Vector{DT}
    q₁::Vector{DT}
    v₁::Vector{DT}

    function InitialGuessODE(interp, rk4, sol, v, Δt, d, periodicity)
        new(interp, rk4, sol, v, Δt, 0, 0, periodicity,
            zeros(DT,d), zeros(DT,d),
            zeros(DT,d), zeros(DT,d))
    end
end

function InitialGuessODE{DT,TT,VT,N}(interp, equation::ODE{DT,TT,VT,N}, Δt::TT; periodicity=[])
    if N > 1
        equ = similar(equation, equation.q₀[:,1])
    else
        equ = equation
    end

    if !(length(periodicity) == equation.d)
        periodicity = zeros(DT, equation.d)
    end

    int = interp(zero(DT), one(DT), Δt, equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuessODE{DT, TT, VT, interp}(int, rk4, sol, equ.v, Δt, equ.d, periodicity)
end

function initialize!{DT,TT,VT,IT}(ig::InitialGuessODE{DT,TT,VT,IT}, t₁::TT, x₁::Vector{DT})
    ig.t₁ = t₁
    set_initial_conditions!(ig.sol, t₁, x₁)
    integrate!(ig.rk4, ig.sol)
    get_data!(ig.sol.q, ig.q₁, 1, 1)
    ig.v(t₁-ig.Δt, ig.q₁, ig.v₁)
end

function update!{DT,TT,VT,IT}(ig::InitialGuessODE{DT,TT,VT,IT}, t₁::TT, x₁::Vector{DT})
    local Δq::DT

    ig.t₀ = ig.t₁
    ig.t₁ = t₁

    simd_copy!(ig.q₁, ig.q₀)
    simd_copy!(ig.v₁, ig.v₀)
    simd_copy!(x₁, ig.q₁)

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

    # compute vector field
    ig.v(t₁, ig.q₁, ig.v₁)
end

function CommonFunctions.evaluate!{DT,TT,VT,IT}(ig::InitialGuessODE{DT,TT,VT,IT}, guess::Vector{DT}, c::TT=one(TT))
    evaluate!(ig.int, ig.q₀, ig.q₁, ig.v₀, ig.v₁, one(TT)+c, guess)
end

function CommonFunctions.evaluate!{DT,TT,VT,IT}(ig::InitialGuessODE{DT,TT,VT,IT},
           guess_q::Vector{DT}, guess_v::Vector{DT}, c::TT=one(TT))

    @assert length(guess_q) == length(guess_v)

    evaluate!(ig.int, ig.q₀, ig.q₁, ig.v₀, ig.v₁, one(TT)+c, guess_q, guess_v)

    if ig.q₀ == ig.q₁
        guess_v .= 0
    end

    nothing
end

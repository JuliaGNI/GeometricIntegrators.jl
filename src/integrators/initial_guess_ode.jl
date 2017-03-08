
using ..CommonFunctions
using ..Interpolation
using ..Tableaus


type InitialGuessODE{DT, TT, VT, IT <: Interpolator}
    int::IT
    rk4::IntegratorERK{DT,TT,VT}
    sol::SolutionODE{DT,TT,2}
    v::VT
    Δt::TT
    q₀::Vector{DT}
    v₀::Vector{DT}
    q₁::Vector{DT}
    v₁::Vector{DT}
end

function InitialGuessODE{DT,TT,VT,N}(interp, equation::ODE{DT,TT,VT,N}, Δt::TT)
    if N > 1
        equ = similar(equation, equation.q₀[:,1])
    else
        equ = equation
    end

    int = interp(zero(DT), one(DT), Δt, equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuessODE{DT, TT, VT, interp}(int, rk4, sol, equ.v, Δt,
                                        zeros(DT, equ.d), zeros(DT, equ.d),
                                        zeros(DT, equ.d), zeros(DT, equ.d))
end

function initialize!{DT,TT,VT,IT}(ig::InitialGuessODE{DT,TT,VT,IT}, t₁::TT, x₁::Vector{DT})
    set_initial_conditions!(ig.sol, t₁, x₁)
    integrate!(ig.rk4, ig.sol)
    get_data!(ig.sol.q, ig.q₁, 1, 1)
    ig.v(t₁-ig.Δt, ig.q₁, ig.v₁)
end

function update!{DT,TT,VT,IT}(ig::InitialGuessODE{DT,TT,VT,IT}, t₁::TT, x₁::Vector{DT})
    simd_copy!(ig.q₁, ig.q₀)
    simd_copy!(ig.v₁, ig.v₀)
    simd_copy!(x₁, ig.q₁)
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

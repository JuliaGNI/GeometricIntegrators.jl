
using ..CommonFunctions
using ..Interpolation
using ..Tableaus


type InitialGuessODE{DT, TT, FT, IT <: Interpolator}
    int::IT
    rk4::IntegratorERK{DT,TT,FT}
    sol::SolutionODE{DT,TT,2}
    f::FT
    Δt::TT
    x₀::Vector{DT}
    f₀::Vector{DT}
    x₁::Vector{DT}
    f₁::Vector{DT}
end

function InitialGuessODE{DT,TT,FT,N}(interp, equation::ODE{DT,TT,FT,N}, Δt::TT)
    if N > 1
        equ = similar(equation, equation.q₀[:,1])
    else
        equ = equation
    end

    int = interp(zero(DT), one(DT), Δt, equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuessODE{DT, TT, FT, interp}(int, rk4, sol, equ.v, Δt,
                                        zeros(DT, equ.d), zeros(DT, equ.d),
                                        zeros(DT, equ.d), zeros(DT, equ.d))
end

function initialize!{DT,TT,FT,IT}(ig::InitialGuessODE{DT,TT,FT,IT}, t₁::TT, x₁::Vector{DT})
    set_initial_conditions!(ig.sol, t₁, x₁)
    integrate!(ig.rk4, ig.sol)
    get_data!(ig.sol.q, ig.x₁, 1, 1)
    ig.f(t₁-ig.Δt, ig.x₁, ig.f₁)
end

function update!{DT,TT,FT,IT}(ig::InitialGuessODE{DT,TT,FT,IT}, t₁::TT, x₁::Vector{DT})
    simd_copy!(ig.x₁, ig.x₀)
    simd_copy!(ig.f₁, ig.f₀)
    simd_copy!(x₁, ig.x₁)
    ig.f(t₁, ig.x₁, ig.f₁)
end

function CommonFunctions.evaluate!{DT,TT,FT,IT}(ig::InitialGuessODE{DT,TT,FT,IT}, guess::Vector{DT}, c::TT=one(TT))
    evaluate!(ig.int, ig.x₀, ig.x₁, ig.f₀, ig.f₁, one(TT)+c, guess)
end

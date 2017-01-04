
using ..CommonFunctions
using ..Interpolation
using ..Tableaus


type InitialGuess{DT, TT, FT, IT <: Interpolator}
    int::IT
    rk4::IntegratorERK{DT,TT,FT}
    sol::SolutionODE{DT,TT}
    f::FT
    Δt::TT
    x₀::Vector{DT}
    f₀::Vector{DT}
    x₁::Vector{DT}
    f₁::Vector{DT}
end

function InitialGuess{DT,TT,FT}(interp, equ::ODE{DT,TT,FT}, Δt::TT)
    int = interp(zero(DT), one(DT), equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuess{DT, TT, FT, int}(int, rk4, sol, equ.v, Δt,
                                  zeros(DT, equ.d), zeros(DT, equ.d),
                                  zeros(DT, equ.d), zeros(DT, equ.d))
end

function initialize!{DT,TT,FT,IT}(ig::InitialGuess{DT,TT,FT,IT}, t₁::TT, x₁::Vector{DT})
    set_initial_conditions!(ig.sol, t₁, x₁)
    integrate!(ig.rk4, ig.sol)
    get_data!(ig.sol.q, ig.x₁, 1, 1)
    ig.f(t₁-ig.Δt, ig.x₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function update!{DT,TT,FT,IT}(ig::InitialGuess{DT,TT,FT,IT}, t₁::TT, x₁::Vector{DT})
    simd_copy!(ig.x₁, ig.x₀)
    simd_copy!(ig.f₁, ig.f₀)
    simd_copy!(x₁, ig.x₁)
    ig.f(t₁, ig.x₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function CommonFunctions.evaluate!{DT,TT,FT,IT}(ig::InitialGuess{DT,TT,FT,IT}, guess::Vector{DT}, c::TT=one(TT))
    evaluate!(ig.int, ig.x₀, ig.x₁, ig.f₀, ig.f₁, one(TT)+c, guess)
end

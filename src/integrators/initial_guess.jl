
include("../interpolation/interpolation.jl")
include("../interpolation/hermite_interpolation.jl")
include("../tableaus/tableaus_erk.jl")

type InitialGuess{DT, TT, FT, IT <: Interpolator}
    int::IT
    rk4::IntegratorERK{DT,TT,FT}
    sol::SolutionODE{DT,TT}
    f::FT
    Δt::TT
    y₀::Vector{DT}
    f₀::Vector{DT}
    y₁::Vector{DT}
    f₁::Vector{DT}
end

function InitialGuess{DT,TT,FT}(int, equ::ODE{DT,TT,FT}, Δt::TT)
    interp = int(zero(DT), one(DT), equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuess{DT, TT, FT, int}(interp, rk4, sol, equ.v, Δt,
                                  zeros(DT, equ.d), zeros(DT, equ.d),
                                  zeros(DT, equ.d), zeros(DT, equ.d))
end

function initialize!{DT,TT,FT,IT}(ig::InitialGuess{DT,TT,FT,IT}, t₁::TT, y₁::Vector{DT})
    set_initial_conditions!(ig.sol, t₁, y₁)
    integrate!(ig.rk4, ig.sol)
    get_data!(ig.sol.q, ig.y₁, 1, 1)
    ig.f(t₁-ig.Δt, ig.y₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function update!{DT,TT,FT,IT}(ig::InitialGuess{DT,TT,FT,IT}, t₁::TT, y₁::Vector{DT})
    simd_copy!(ig.y₁, ig.y₀)
    simd_copy!(ig.f₁, ig.f₀)
    simd_copy!(y₁, ig.y₁)
    ig.f(t₁, ig.y₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function evaluate{DT,TT,FT,IT}(ig::InitialGuess{DT,TT,FT,IT}, guess::Vector{DT}, c::TT=one(TT))
    evaluate(ig.int, ig.y₀, ig.y₁, ig.f₀, ig.f₁, one(TT)+c, guess)
end

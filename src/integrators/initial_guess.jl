
include("../interpolation/interpolation.jl")
include("../interpolation/hermite_interpolation.jl")
include("../tableaus/tableaus_erk.jl")

type InitialGuess{T, FT, IT <: Interpolator}
    int::IT
    rk4::IntegratorERK{T,FT}
    sol::SolutionODE{T}
    f::FT
    Δt::T
    y₀::Vector{T}
    f₀::Vector{T}
    y₁::Vector{T}
    f₁::Vector{T}
end

function InitialGuess{T,FT}(int, equ::ODE{T,FT}, Δt::T)
    interp = int(zero(T), one(T), equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuess{T, FT, int}(interp, rk4, sol, equ.f, Δt,
                             zeros(T, equ.d), zeros(T, equ.d),
                             zeros(T, equ.d), zeros(T, equ.d))
end

function initialize!{T,FT,IT}(ig::InitialGuess{T,FT,IT}, t₁::T, y₁::Vector{T})
    set_initial_conditions!(ig.sol, t₁, y₁)
    integrate!(ig.rk4, ig.sol)
    simd_copy_xy_first!(ig.y₁, ig.sol, 1, 1)
    ig.f(t₁-ig.Δt, ig.y₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function update!{T,FT,IT}(ig::InitialGuess{T,FT,IT}, t₁::T, y₁::Vector{T})
    simd_copy!(ig.y₁, ig.y₀)
    simd_copy!(ig.f₁, ig.f₀)
    simd_copy!(y₁, ig.y₁)
    ig.f(t₁, ig.y₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function evaluate{T,FT,IT}(ig::InitialGuess{T,FT,IT}, guess::Vector{T}, c::T=one(T))
    evaluate(ig.int, ig.y₀, ig.y₁, ig.f₀, ig.f₁, one(T)+c, guess)
end

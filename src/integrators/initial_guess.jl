
include("../interpolation/interpolation.jl")
include("../interpolation/hermite_interpolation.jl")
include("../tableaus/tableaus_erk.jl")

type InitialGuess{T, I <: Interpolator}
    int::I
    rk4::IntegratorERK{T}
    sol::SolutionODE{T}
    f::Function
    Δt::T
    y₀::Vector{T}
    f₀::Vector{T}
    y₁::Vector{T}
    f₁::Vector{T}
end

function InitialGuess{T}(int, equ::ODE{T}, Δt::T)
    interp = int(zero(T), one(T), equ.d)
    rk4 = IntegratorERK(equ, getTableauERK4(), -Δt)
    sol = SolutionODE(equ, Δt, 1)
    InitialGuess{T, int}(interp, rk4, sol, equ.f, Δt,
                         zeros(T, equ.d), zeros(T, equ.d),
                         zeros(T, equ.d), zeros(T, equ.d))
end

function initialize!{T}(ig::InitialGuess{T}, t₁::T, y₁::Vector{T})
    set_initial_conditions!(ig.sol, t₁, y₁)
    integrate!(ig.rk4, ig.sol)
    simd_copy_xy_first!(ig.y₁, ig.sol, 1, 1)
    ig.f(t₁-ig.Δt, ig.y₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function update!{T}(ig::InitialGuess{T}, t₁::T, y₁::Vector{T})
    simd_copy!(ig.y₁, ig.y₀)
    simd_copy!(ig.f₁, ig.f₀)
    simd_copy!(y₁, ig.y₁)
    ig.f(t₁, ig.y₁, ig.f₁)
    simd_scale!(ig.f₁, ig.Δt)
end

function evaluate{T}(ig::InitialGuess{T}, guess::Vector{T}, c::T=one(T))
    evaluate(ig.int, ig.y₀, ig.y₁, ig.f₀, ig.f₁, one(T)+c, guess)
end

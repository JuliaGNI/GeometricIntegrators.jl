
using ..CommonFunctions

immutable HermiteInterpolation{T} <: Interpolator{T}

    x₀::T
    x₁::T
    a₀::T
    a₁::T
    b₀::T
    b₁::T

    num::Vector{T}

    function HermiteInterpolation(x₀, x₁, d)
        a₀ = 2/(x₀-x₁)
        a₁ = 2/(x₁-x₀)

        b₀ = one(T)/(x₀-x₁)^2
        b₁ = one(T)/(x₁-x₀)^2

        new(x₀, x₁, a₀, a₁, b₀, b₁, zeros(T, d))
    end
end

function HermiteInterpolation{T}(x₀::T, x₁::T, d::Int)
    HermiteInterpolation{T}(x₀, x₁, d)
end


function CommonFunctions.evaluate!{T}(int::HermiteInterpolation{T}, y₀::Vector{T}, y₁::Vector{T}, f₀::Vector{T}, f₁::Vector{T}, x::T, y::Vector{T})
    local d₀::T
    local d₁::T
    local c₀::T
    local c₁::T
    local den::T

    # Interpolate y values at required locations
    if x == int.x₀
        simd_copy!(y₀, y)
    elseif x == int.x₁
        simd_copy!(y₁, y)
    else
        d₀ = one(T)/(x-int.x₀)
        c₁ = d₀*int.b₀
        c₀ = c₁*(d₀-int.a₀)

        simd_copy_scale!(c₀, y₀, int.num)
        simd_axpy!(c₁, f₀, int.num)
        den = c₀

        d₁ = one(T)/(x-int.x₁)
        c₁ = d₁*int.b₁
        c₀ = c₁*(d₁-int.a₁)

        simd_axpy!(c₀, y₁, int.num)
        simd_axpy!(c₁, f₁, int.num)
        den += c₀

        simd_copy_scale!(one(T)/den, int.num, y)
    end
end

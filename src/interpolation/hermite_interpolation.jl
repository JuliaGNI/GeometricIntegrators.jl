@doc raw"""
# Hermite's Interpolating Polynomials

Here, we implement a two point Hermite interpolation function which passes
through the function and its first derivative for the interval ``[0,1]``.
The polynomial is determined by four constraint equations, matching the
function and its derivative at the points ``0`` and ``1``.

These inter- and extrapolation routines are exclusively used for computing
initial guesses and are usually not called directly by the user.

Start by defining the 3rd degree polynomial and its derivative by
```math
\begin{aligned}
g(x) &= a_0 + a_1 x + a_2 x^2 + a_3 x^3 , \\
g'(x) &= a_1 + 2 a_2 x + 3 a_3 x^2 ,
\end{aligned}
```
and apply the constraints
```math
\begin{aligned}
g(0) &= f_0 & & \Rightarrow & a_0 &= f_0 , \\
g(1) &= f_1 & & \Rightarrow & a_0 + a_1 + a_2 + a_3 &= f_1 , \\
g'(0) &= f'_0 & & \Rightarrow & a_1 &= f'_0 , \\
g'(1) &= f'_1 & & \Rightarrow & a_1 + 2 a_2 + 3 a_3 &= f'_1 . \\
\end{aligned}
```
Solving for ``a_0, a_1, a_2, a_3`` leads to
```math
\begin{aligned}
a_0 &= f_0 , &
a_1 &= f'_0 , &
a_2 &= - 3 f_0 + 3 f_1 - 2 f'_0 - f'_1 , &
a_3 &= 2 f_0 - 2 f_1 + f'_0 + f'_1 ,
\end{aligned}
```
so that the polynomial ``g(x)`` reads
```math
g(x) = f_0 + f'_0 x + (- 3 f_0 + 3 f_1 - 2 f'_0 - f'_1) x^2 + (2 f_0 - 2 f_1 + f'_0 + f'_1) x^3 .
```
The function and derivative values can be factored out, so that ``g(x)`` can be rewritten as
```math
g(x) = f_0 (1 - 3 x^2 + 2 x^3) + f_1 (3 x^2 - 2 x^3) + f'_0 (x - 2 x^2 + x^3) + f'_1 (- x^2 + x^3) ,
```
or in generic form as
```math
g(x) = f_0 a_0(x) + f_1 a_1(x) + f'_0 b_0(x) + f'_1 b_1(x) ,
```
with basis functions
```math
\begin{aligned}
a_0 (x) &= 1 - 3 x^2 + 2 x^3 , &
b_0 (x) &= x - 2 x^2 + x^3 , \\
a_1 (x) &= 3 x^2 - 2 x^3 , &
b_1 (x) &= - x^2 + x^3 .
\end{aligned}
```
The derivative ``g'(x)`` accordingly reads
```math
g'(x) = f_0 a'_0(x) + f_1 a'_1(x) + f'_0 b'_0(x) + f'_1 b'_1(x) ,
```
with
```math
\begin{aligned}
a'_0 (x) &= - 6 x + 6 x^2 , &
b'_0 (x) &= 1 - 4 x + 3 x^2 , \\
a'_1 (x) &= 6 x - 6 x^2 , &
b'_1 (x) &= - 2 x + 3 x^2 .
\end{aligned}
```
The basis functions ``a_0``and ``a_1`` are associated with the function
values at ``x_0`` and ``x_1``, respectively, while the basis functions
``b_0`` and ``b_1`` are associated with the derivative values at
``x_0`` and ``x_1``.
The basis functions satisfy the following relations,
```math
\begin{aligned}
a_i (x_j) &= \delta_{ij} , &
b_i (x_j) &= 0 , &
a'_i (x_j) &= 0 , &
b'_i (x_j) &= \delta_{ij} , &
i,j &= 0, 1 ,
\end{aligned}
```
where ``\delta_{ij}`` denotes the Kronecker-delta, so that
```math
\begin{aligned}
g(0) &= f_0 , &
g(1) &= f_1 , &
g'(0) &= f'_0 , &
g'(1) &= f'_1 .
\end{aligned}
```
"""
struct HermiteInterpolation{T} <: Interpolator{T}
    x₀::T
    x₁::T
    Δx::T

    function HermiteInterpolation{T}(x₀, x₁, Δx) where {T}
        new(x₀, x₁, Δx)
    end
end

function HermiteInterpolation(x₀::T, x₁::T, Δx::T) where {T}
    HermiteInterpolation{T}(x₀, x₁, Δx)
end


function Common.evaluate!(int::HermiteInterpolation{TT}, y₀::AbstractArray{DT}, y₁::AbstractArray{DT}, f₀::AbstractArray{DT}, f₁::AbstractArray{DT}, x::TT, y::AbstractArray{DT}) where {DT,TT}
    local a₀::TT
    local a₁::TT
    local b₀::TT
    local b₁::TT

    # Interpolate y values at required locations
    if x == int.x₀
        y .= y₀
    elseif x == int.x₁
        y .= y₁
    else
        a₁ = 3x^2 - 2x^3
        a₀ = 1 - a₁
        b₁ = x^2*(x-1)
        b₀ = x*(1-x)+b₁
        y .= a₀ .* y₀ .+ a₁ .* y₁ .+ b₀ .* int.Δx .* f₀ .+ b₁ .* int.Δx .* f₁
    end
end

function Common.evaluate!(int::HermiteInterpolation{TT}, y₀::AbstractArray{DT}, y₁::AbstractArray{DT}, f₀::AbstractArray{DT}, f₁::AbstractArray{DT}, x::TT, y::AbstractArray{DT}, f::AbstractArray{DT}) where {DT,TT}
    local a₀::TT
    local a₁::TT
    local b₀::TT
    local b₁::TT

    evaluate!(int, y₀, y₁, f₀, f₁, x, y)

    # Interpolate f values at required locations
    if x == int.x₀
        f .= f₀
    elseif x == int.x₁
        f .= f₁
    else
        a₁ = (6x - 6x^2) / int.Δx
        a₀ = - a₁
        b₁ = x*(3x-2)
        b₀ = 1-2x+b₁
        f .= a₀ .* y₀ .+ a₁ .* y₁ .+ b₀ .* f₀ .+ b₁ .* f₁
    end
end

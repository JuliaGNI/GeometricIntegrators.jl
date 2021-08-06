@doc raw"""
# Hermite's Interpolating Polynomials

Implements a two point Hermite inter-/extrapolation function which passes
through the function and its first derivative for the interval ``[0,1]``.
The polynomial is determined by four constraint equations, matching the
function and its derivative at the points ``0`` and ``1``.

Call with one of the following methods
```julia
_hermite_extrapolation!(t₀, t₁, x₀, x₁, ẋ₀, ẋ₁, t, x)
_hermite_extrapolation!(t₀, t₁, x₀, x₁, ẋ₀, ẋ₁, t, x, ẋ)
_hermite_extrapolation!(v, t₀, t₁, x₀, x₁, t, x)
_hermite_extrapolation!(v, t₀, t₁, x₀, x₁, t, x, ẋ)
```

where

* `t₀`: first  sample time $t_0$
* `t₁`: second sample time $t_1$
* `x₀`: first  solution value $x_0 = x(t_0)$
* `x₁`: second solution value $x_1 = x(t_1)$
* `ẋ₀`: first  vector field value $ẋ_0 = v(t_0, x(t_0))$
* `ẋ₁`: second vector field value $ẋ_1 = v(t_1, x(t_1))$
* `v`:  function to compute vector field with signature `v(t,x,ẋ)`
* `t`:  time $t$ to extrapolate
* `x`:  extrapolated solution value $x(t)$
* `v`:  extrapolated vector field value $ẋ(t)$


#### Derivation

The interpolation works as follows:
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
function _hermite_extrapolation! end


function _hermite_extrapolation!(t₀::TT, t₁::TT, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}, ẋ₀::AbstractArray{DT}, ẋ₁::AbstractArray{DT}, t::TT, x::AbstractArray{DT}) where {DT,TT}
    local a₀::TT
    local a₁::TT
    local b₀::TT
    local b₁::TT
    local Δt = t₁ - t₀
    local s = (t - t₀) / Δt

    # Interpolate x at t
    if t == t₀
        x .= x₀
    elseif t == t₁
        x .= x₁
    else
        a₁ = 3s^2 - 2s^3
        a₀ = 1 - a₁
        b₁ = s^2*(s-1)
        b₀ = s*(1-s)+b₁
        x .= a₀ .* x₀ .+ a₁ .* x₁ .+ b₀ .* Δt .* ẋ₀ .+ b₁ .* Δt .* ẋ₁
    end

    return x
end

function _hermite_extrapolation!(t₀::TT, t₁::TT, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}, ẋ₀::AbstractArray{DT}, ẋ₁::AbstractArray{DT}, t::TT, x::AbstractArray{DT}, ẋ::AbstractArray{DT}) where {DT,TT}
    local a₀::TT
    local a₁::TT
    local b₀::TT
    local b₁::TT
    local Δt = t₁ - t₀
    local s = (t - t₀) / Δt

    _hermite_extrapolation!(t₀, t₁, x₀, x₁, ẋ₀, ẋ₁, t, x)

    # Interpolate ẋ at t
    if t == t₀
        ẋ .= ẋ₀
    elseif t == t₁
        ẋ .= ẋ₁
    else
        a₁ = (6s - 6s^2) / Δt
        a₀ = - a₁
        b₁ = s*(3s-2)
        b₀ = 1-2s+b₁
        ẋ .= a₀ .* x₀ .+ a₁ .* x₁ .+ b₀ .* ẋ₀ .+ b₁ .* ẋ₁
    end

    return (x, ẋ)
end

function _get_velocities(v::Function, t₀::TT, t₁::TT, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}) where {DT,TT}
    ẋ₀ = zero(x₀)
    ẋ₁ = zero(x₁)
    v(t₀, x₀, ẋ₀)
    v(t₁, x₁, ẋ₁)
    return (ẋ₀, ẋ₁)
end

function _hermite_extrapolation!(v::Function, t₀::TT, t₁::TT, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}, t::TT, x::AbstractArray{DT}) where {DT,TT}
    _hermite_extrapolation!(t₀, t₁, x₀, x₁, _get_velocities(v, t₀, t₁, x₀, x₁)..., t, x)
end

function _hermite_extrapolation!(v::Function, t₀::TT, t₁::TT, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}, t::TT, x::AbstractArray{DT}, ẋ::AbstractArray{DT}) where {DT,TT}
    _hermite_extrapolation!(t₀, t₁, x₀, x₁, _get_velocities(v, t₀, t₁, x₀, x₁)..., t, x, ẋ)
end



struct HermiteExtrapolation{T} <: Extrapolation
    t₀::T
    t₁::T
    Δt::T

    function HermiteExtrapolation{T}(t₀, t₁) where {T}
        new(t₀, t₁, t₁-t₀)
    end
end

function HermiteExtrapolation(t₀::T, t₁::T) where {T}
    HermiteExtrapolation{T}(t₀, t₁)
end


function GeometricBase.evaluate!(int::HermiteExtrapolation{TT}, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}, ẋ₀::AbstractArray{DT}, ẋ₁::AbstractArray{DT}, t::TT, x::AbstractArray{DT}) where {DT,TT}
    _hermite_extrapolation!(int.t₀, int.t₁, x₀, x₁, ẋ₀, ẋ₁, t, x)
end

function GeometricBase.evaluate!(int::HermiteExtrapolation{TT}, x₀::AbstractArray{DT}, x₁::AbstractArray{DT}, ẋ₀::AbstractArray{DT}, ẋ₁::AbstractArray{DT}, t::TT, x::AbstractArray{DT}, ẋ::AbstractArray{DT}) where {DT,TT}
    _hermite_extrapolation!(int.t₀, int.t₁, x₀, x₁, ẋ₀, ẋ₁, t, x, ẋ)
end

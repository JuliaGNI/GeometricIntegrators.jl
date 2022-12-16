@doc raw"""
# Hermite's Interpolating Polynomials

Implements a two point Hermite inter-/extrapolation function which passes
through the function and its first derivative for the interval ``[0,1]``.
The polynomial is determined by four constraint equations, matching the
function and its derivative at the points ``0`` and ``1``.

Call with one of the following methods
```julia
extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, t, x, HermiteExtrapolation())
extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, t, x, ẋ, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, v, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, ẋ, v, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, problem, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, ẋ, problem, HermiteExtrapolation())
```

where

* `t₀`: first  sample time $t_0$
* `x₀`: first  solution value $x_0 = x(t_0)$
* `ẋ₀`: first  vector field value $ẋ_0 = v(t_0, x(t_0))$
* `t₁`: second sample time $t_1$
* `x₁`: second solution value $x_1 = x(t_1)$
* `ẋ₁`: second vector field value $ẋ_1 = v(t_1, x(t_1))$
* `t`:  time $t$ to extrapolate
* `x`:  extrapolated solution value $x(t)$
* `ẋ`:  extrapolated vector field value $ẋ(t)$
* `v`:  function to compute vector field with signature `v(ẋ,t,x)`
* `problem`: [`ODEProblem`](@ref) whose vector field to use


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
struct HermiteExtrapolation end


function extrapolate!(t₀::TT, x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
                      tᵢ::TT, xᵢ::AbstractArray{DT},
                      ::HermiteExtrapolation) where {DT,TT}
    local a₀::TT
    local a₁::TT
    local b₀::TT
    local b₁::TT
    local Δt = t₁ - t₀
    local s = (tᵢ - t₀) / Δt

    # Interpolate x at t
    if tᵢ == t₀
        xᵢ .= x₀
    elseif tᵢ == t₁
        xᵢ .= x₁
    else
        a₁ = 3s^2 - 2s^3
        a₀ = 1 - a₁
        b₁ = s^2*(s-1)
        b₀ = s*(1-s)+b₁
        xᵢ .= a₀ .* x₀ .+ a₁ .* x₁ .+ b₀ .* Δt .* ẋ₀ .+ b₁ .* Δt .* ẋ₁
    end

    return xᵢ
end

function extrapolate!(t₀::TT, x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
                      tᵢ::TT, xᵢ::AbstractArray{DT}, ẋᵢ::AbstractArray{DT},
                      extrap::HermiteExtrapolation) where {DT,TT}
    local a₀::TT
    local a₁::TT
    local b₀::TT
    local b₁::TT
    local Δt = t₁ - t₀
    local s = (tᵢ - t₀) / Δt

    extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, tᵢ, xᵢ, extrap)

    # Interpolate ẋ at t
    if tᵢ == t₀
        ẋᵢ .= ẋ₀
    elseif tᵢ == t₁
        ẋᵢ .= ẋ₁
    else
        a₁ = (6s - 6s^2) / Δt
        a₀ = - a₁
        b₁ = s*(3s-2)
        b₀ = 1-2s+b₁
        ẋᵢ .= a₀ .* x₀ .+ a₁ .* x₁ .+ b₀ .* ẋ₀ .+ b₁ .* ẋ₁
    end

    return (xᵢ, ẋᵢ)
end


function _vectorfield(t::TT, x::AbstractArray{DT}, v::Callable) where {DT,TT}
    ẋ = zero(x)
    v(ẋ, t, x)
    return ẋ
end

function extrapolate!(t₀::TT, x₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT},
                      tᵢ::TT, xᵢ::AbstractArray{DT},
                      v::Function, extrap::HermiteExtrapolation) where {DT,TT}
    ẋ₀ = _vectorfield(t₀, x₀, v)
    ẋ₁ = _vectorfield(t₁, x₁, v)
    extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, tᵢ, xᵢ, extrap)
end

function extrapolate!(t₀::TT, x₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT},
                      tᵢ::TT, xᵢ::AbstractArray{DT}, ẋᵢ::AbstractArray{DT}, 
                      v::Function, extrap::HermiteExtrapolation) where {DT,TT}
    ẋ₀ = _vectorfield(t₀, x₀, v)
    ẋ₁ = _vectorfield(t₁, x₁, v)
    extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, tᵢ, xᵢ, ẋᵢ, extrap)
end

function extrapolate!(t₀::TT, x₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT},
                      tᵢ::TT, xᵢ::AbstractArray{DT},
                      problem::ODEProblem,
                      extrap::HermiteExtrapolation) where {DT,TT}
    extrapolate!(t₀, x₀, t₁, x₁, tᵢ, xᵢ, functions(problem).v, extrap)
end

function extrapolate!(t₀::TT, x₀::AbstractArray{DT},
                      t₁::TT, x₁::AbstractArray{DT},
                      tᵢ::TT, xᵢ::AbstractArray{DT}, ẋᵢ::AbstractArray{DT}, 
                      problem::ODEProblem,
                      extrap::HermiteExtrapolation) where {DT,TT}
    extrapolate!(t₀, x₀, t₁, x₁, tᵢ, xᵢ, ẋᵢ, functions(problem).v, extrap)
end

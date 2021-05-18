
@doc raw"""
Euler extrapolation method with arbitrary order p.

Solves the ordinary differential equation
```math
\begin{aligned}
\dot{x} &= v(t, x) , &
x(t_0) &= x_0 ,
\end{aligned}
```
for $x_1 = x(t_1)$.

Call with
```julia
_euler_extrapolation_ode!(v, t₀, t₁, x₀, x₁, s)
```

where

* `v`:  function to compute vector field with signature `v(t,x,ẋ)`
* `t₀`: initial time
* `t₁`: final   time
* `x₀`: initial value $x_0 = x(t_0)$
* `x₁`: final   value $x_1 = x(t_1)$
* `s`:  number of interpolations (order $p=s+1$)

"""
function _euler_extrapolation_ode!(v::Function, t₀::TT, t₁::TT, x₀::AbstractVector{DT}, x₁::AbstractVector{DT}, s::Int) where {DT,TT}
    @assert size(x₀) == size(x₁)

    F   = collect(1:(s+1))
    Δt  = t₁ - t₀
    σ   = Δt ./ F
    pts = repeat(x₀, outer = [1, s+1])

    local xᵢ = zero(x₀)
    local vᵢ = zero(x₀)

    for i in F
        for _ in 1:(F[i]-1)
            tᵢ = t₀ + σ[i]
            for k in axes(pts,1)
                xᵢ[k] = pts[k,i]
            end
            v(tᵢ, xᵢ, vᵢ)
            for k in axes(pts,1)
                pts[k,i] += σ[i] * vᵢ[k]
            end
        end
    end

    aitken_neville!(σ, pts, zero(TT), x₁)
    return x₁
end


# struct EulerExtrapolation <: Extrapolation
#     s::Int
# end

# function Common.evaluate!(extrap::EulerExtrapolation, v::Function, t₀::TT, t₁::TT, x₀::AbstractVector{DT}, x₁::AbstractVector{DT}) where {DT,TT}
#     _euler_extrapolation_ode!(v, t₀, t₁, x₀, x₁, extrap.s)
# end

# function Common.evaluate!(extrap::EulerExtrapolation, v::Function, f::Function, t₀::TT, t₁::TT, q₀::AbstractVector{DT}, q₁::AbstractVector{DT}, p₀::AbstractVector{DT}, p₁::AbstractVector{DT}) where {DT,TT}
#     _euler_extrapolation_pode!(v, t₀, t₁, q₀, q₁, p₀, p₁, extrap.s)
# end


struct EulerExtrapolationODE{VT} <: Extrapolation
    s::Int
    v::VT

    function EulerExtrapolationODE{VT}(v, s) where {VT}
        new(s,v)
    end
end

function EulerExtrapolationODE(v::VT, s::Int) where {VT}
    EulerExtrapolationODE{VT}(v, s)
end

function EulerExtrapolation(equ::ODE, s::Int)
    EulerExtrapolationODE(_get_v(equ), s)
end

function Common.evaluate!(extrap::EulerExtrapolationODE, t₀::TT, t₁::TT, x₀::AbstractVector{DT}, x₁::AbstractVector{DT}) where {DT,TT}
    _euler_extrapolation_ode!(extrap.v, t₀, t₁, x₀, x₁, extrap.s)
end


# struct EulerExtrapolationIODE{VT,FT} <: Extrapolation
#     s::Int
#     v::VT
#     f::FT
    
#     function EulerExtrapolationIODE{VT,FT}(v, f, s) where {VT,FT}
#         new(s,v,f)
#     end
# end

# function EulerExtrapolationIODE(v::VT, f::FT, s::Int) where {VT,FT}
#     EulerExtrapolationIODE{VT,FT}(v, f, s)
# end

# function EulerExtrapolation(equ::IODE, s::Int)
#     EulerExtrapolationIODE(_get_v(equ), _get_f(equ), s)
# end

# function Common.evaluate!(extrap::EulerExtrapolationIODE, t₀::TT, t₁::TT, q₀::AbstractVector{DT}, q₁::AbstractVector{DT}, p₀::AbstractVector{DT}, p₁::AbstractVector{DT}) where {DT,TT}
#     _euler_extrapolation_iode!(extrap.v, t₀, t₁, q₀, q₁, p₀, p₁, extrap.s)
# end


# struct EulerExtrapolationPODE{VT,FT} <: Extrapolation
#     s::Int
#     v::VT
#     f::FT
    
#     function EulerExtrapolationPODE{VT,FT}(v, f, s) where {VT,FT}
#         new(s,v,f)
#     end
# end

# function EulerExtrapolationPODE(v::VT, f::FT, s::Int) where {VT,FT}
#     EulerExtrapolationPODE{VT,FT}(v, f, s)
# end

# function EulerExtrapolation(equ::PODE, s::Int)
#     EulerExtrapolationPODE(_get_v(equ), _get_f(equ), s)
# end

# function Common.evaluate!(extrap::EulerExtrapolationPODE, t₀::TT, t₁::TT, q₀::AbstractVector{DT}, q₁::AbstractVector{DT}, p₀::AbstractVector{DT}, p₁::AbstractVector{DT}) where {DT,TT}
#     _euler_extrapolation_pode!(extrap.v, t₀, t₁, q₀, q₁, p₀, p₁, extrap.s)
# end

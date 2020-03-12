@doc raw"""
PathIntegralLinear is a path integral along a linear path

```math
\phi (\tau; q^-, q^+) = (1-\tau) q^- + \tau q^+ .
```
"""
struct PathIntegralLinear <: PathIntegral end


# function evaluate(path::PathIntegralLinear{T}, i::Int, τ::T) where {T}
#     if i == 1
#         return one(T)-τ
#     elseif i == 2
#         return τ
#     else
#         error("PathIntegralLinear has only two basis functions, not ", k, ".")
#     end
# end

# function derivative(path::PathIntegralLinear{T}, i::Int, τ::T) where {T}
#     if i == 1
#         return -one(T)
#     elseif i == 2
#         return +one(T)
#     else
#         error("PathIntegralLinear has only two basis functions, not ", k, ".")
#     end
# end

evaluate_l(path::PathIntegralLinear, τ::T) where {T} = one(T)-τ
evaluate_r(path::PathIntegralLinear, τ::T) where {T} = τ

derivative_l(path::PathIntegralLinear, τ::T) where {T} = -one(T)
derivative_r(path::PathIntegralLinear, τ::T) where {T} = +one(T)

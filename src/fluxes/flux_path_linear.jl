"""
FluxPathLinear is a linear path

```math
\\phi (\\tau; q^-, q^+) = (1-\\tau) q^- + \\tau q^+ .
```
"""
struct FluxPathLinear <: FluxPath end


# function evaluate(path::FluxPathLinear{T}, i::Int, τ::T) where {T}
#     if i == 1
#         return one(T)-τ
#     elseif i == 2
#         return τ
#     else
#         error("FluxPathLinear has only two basis functions, not ", k, ".")
#     end
# end

# function derivative(path::FluxPathLinear{T}, i::Int, τ::T) where {T}
#     if i == 1
#         return -one(T)
#     elseif i == 2
#         return +one(T)
#     else
#         error("FluxPathLinear has only two basis functions, not ", k, ".")
#     end
# end

evaluate_l(path::FluxPathLinear, τ::T) where {T} = one(T)-τ
evaluate_r(path::FluxPathLinear, τ::T) where {T} = τ

derivative_l(path::FluxPathLinear, τ::T) where {T} = -one(T)
derivative_r(path::FluxPathLinear, τ::T) where {T} = +one(T)

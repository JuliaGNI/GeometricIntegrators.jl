"""
FluxPathTrigonometric is a cos^2/sin^2 path

```math
\\phi (\\tau; q^-, q^+) = \\cos^2 (\\pi \\tau / 2) q^- + \\sin^2 (\\pi \\tau / 2) q^+ .
```
"""
struct FluxPathTrigonometric <: FluxPath end

evaluate_l(path::FluxPathTrigonometric, τ::T) where {T} = cos(π/2*τ)^2
evaluate_r(path::FluxPathTrigonometric, τ::T) where {T} = sin(π/2*τ)^2

derivative_l(path::FluxPathTrigonometric, τ::T) where {T} = -π*sin(π/2*τ)*cos(π/2*τ)
derivative_r(path::FluxPathTrigonometric, τ::T) where {T} = +π*sin(π/2*τ)*cos(π/2*τ)

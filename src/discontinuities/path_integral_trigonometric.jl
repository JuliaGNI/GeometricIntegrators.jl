@doc raw"""
PathIntegralTrigonometric is a path integral along a cos^2/sin^2 path

```math
\phi (\tau; q^-, q^+) = \cos^2 (\pi \tau / 2) q^- + \sin^2 (\pi \tau / 2) q^+ .
```
"""
struct PathIntegralTrigonometric <: PathIntegral end

evaluate_l(path::PathIntegralTrigonometric, τ::T) where {T} = cos(π/2*τ)^2
evaluate_r(path::PathIntegralTrigonometric, τ::T) where {T} = sin(π/2*τ)^2

derivative_l(path::PathIntegralTrigonometric, τ::T) where {T} = -π*sin(π/2*τ)*cos(π/2*τ)
derivative_r(path::PathIntegralTrigonometric, τ::T) where {T} = +π*sin(π/2*τ)*cos(π/2*τ)

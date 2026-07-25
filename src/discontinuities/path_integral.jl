@doc raw"""
Abstract supertype of the paths used to integrate a Lagrangian across a jump of a
discontinuous discrete trajectory.

A concrete `PathIntegral` provides the two functions ``\varphi^{-}`` and
``\varphi^{+}`` that build the connecting path from the one-sided limits,
```math
\Phi (\tau; q^-, q^+) = q^- \, \varphi^{-} (\tau) + q^+ \, \varphi^{+} (\tau) ,
\qquad
\varphi^{-} (0) = \varphi^{+} (1) = 1 ,
\qquad
\varphi^{-} (1) = \varphi^{+} (0) = 0 ,
```
through the interface

```julia
evaluate_l(path, τ)     # φ⁻(τ)
evaluate_r(path, τ)     # φ⁺(τ)
derivative_l(path, τ)   # dφ⁻/dτ (τ)
derivative_r(path, τ)   # dφ⁺/dτ (τ)
```

Implementations: [`PathIntegralLinear`](@ref), [`PathIntegralTrigonometric`](@ref).
Paired with a quadrature rule in a [`Discontinuity`](@ref).
"""
abstract type PathIntegral end

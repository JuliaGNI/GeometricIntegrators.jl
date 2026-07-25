@doc raw"""
A discontinuity of the discrete trajectory, together with the means of integrating a
Lagrangian across it.

Discontinuous Galerkin variational integrators approximate the trajectory by a
piecewise polynomial that need not be continuous at the interval boundaries. The
resulting non-conservative product of the one-form and the jump is given a meaning by
integrating along a path ``\Phi`` connecting the two one-sided limits ``q^-`` and
``q^+``,
```math
\int \limits_0^1 \vartheta \big( \Phi(\tau; q^-, q^+) \big) \cdot \frac{d \Phi}{d\tau} (\tau; q^-, q^+) \, d\tau
\approx \sum \limits_{i=1}^{\sigma} \beta_i \, \vartheta \big( \Phi (\gamma_i) \big) \cdot \Phi' (\gamma_i) ,
```
so a `Discontinuity` pairs

* `path`: the connecting path, a [`PathIntegral`](@ref) such as
  [`PathIntegralLinear`](@ref) or [`PathIntegralTrigonometric`](@ref), and
* `quadrature`: the quadrature rule with ``\sigma`` nodes ``\gamma_i`` and weights
  ``\beta_i`` used to evaluate that integral.

Used by `DGVIPI`.

### Constructor

```
Discontinuity(path::PathIntegral, quadrature::QuadratureRule)
```

### Example

```julia
Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))
```

which integrates along a straight line between the one-sided limits with the
two-point Lobatto rule, i.e. a trapezoidal flux.
"""
struct Discontinuity{T<:AbstractFloat, PT<:PathIntegral, QN}
    path::PT
    quadrature::QuadratureRule{T,QN}
end

function Discontinuity(path::PT, quadrature::QuadratureRule{T,QN}) where {T,PT,QN}
    Discontinuity{T,PT,QN}(path, quadrature)
end

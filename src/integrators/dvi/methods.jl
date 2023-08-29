
abstract type DVIMethod <: LODEMethod end


"Symplectic Euler-A Degenerate Variational Integrator."
struct DVIA <: DVIMethod end

"Symplectic Euler-B Degenerate Variational Integrator."
struct DVIB <: DVIMethod end

"Midpoint Degenerate Variational Integrator."
struct CMDVI <: DVIMethod end

"Trapezoidal Degenerate Variational Integrator."
struct CTDVI <: DVIMethod end

GeometricBase.order(::Union{DVIA, Type{DVIA}}) = 1
GeometricBase.order(::Union{DVIB, Type{DVIB}}) = 1
GeometricBase.order(::Union{CMDVI, Type{CMDVI}}) = 2
GeometricBase.order(::Union{CTDVI, Type{CTDVI}}) = 2

issymmetric(::Union{DVIA, Type{<:DVIA}}) = false
issymmetric(::Union{DVIB, Type{<:DVIB}}) = false
issymmetric(::Union{CMDVI, Type{<:CMDVI}}) = true
issymmetric(::Union{CTDVI, Type{<:CTDVI}}) = true

issymplectic(::Union{DVIMethod, Type{<:DVIMethod}}) = true
isexplicit(::Union{DVIMethod, Type{<:DVIMethod}}) = false
isimplicit(::Union{DVIMethod, Type{<:DVIMethod}}) = true


"""
Degenerate Variational Runge-Kutta Methods

```
DVRK(tableau::Tableau)
DVRK(method::RKMethod)
```
"""
struct DVRK{TT} <: DVIMethod
    tableau::TT

    function DVRK(tableau::TT) where {TT <: Tableau}
        new{TT}(tableau)
    end
end

DVRK(method::RKMethod, args...; kwargs...) = DVRK(tableau(method))

GeometricBase.tableau(method::DVRK) = method.tableau
GeometricBase.order(method::DVRK) = order(tableaus(method))
isexplicit(method::DVRK) = false
isimplicit(method::DVRK) = true
issymmetric(method::DVRK) = issymmetric(tableaus(method))
issymplectic(method::DVRK) = issymplectic(tableaus(method))

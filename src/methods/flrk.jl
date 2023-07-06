"""
Formal Lagrangian Runge-Kutta Method

```
FLRK(tableau::Tableau)
FLRK(method::RKMethod)
```
"""
struct FLRK{TT} <: RKMethod
    tableau::TT

    function FLRK(tableau::TT) where {TT <: Tableau}
        new{TT}(tableau)
    end
end

FLRK(method::RKMethod, args...; kwargs...) = FLRK(tableau(method))

GeometricBase.tableau(method::FLRK) = method.tableau

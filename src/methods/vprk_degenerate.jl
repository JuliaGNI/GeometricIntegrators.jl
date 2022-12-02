"""
Variational Partitioned Runge-Kutta Method for Degenerate Lagrangians

```
DegenerateVPRK(method::VPRK)
DegenerateVPRK(tableau::AbstractTableau, args...; kwargs...)
```
"""
struct DegenerateVPRK{VT} <: VPRKMethod
    vprk::VT

    function DegenerateVPRK(vprk::VT) where {VT <: VPRKMethod}
        new{VT}(vprk)
    end
end

DegenerateVPRK(tableau::AbstractTableau, args...; kwargs...) = DegenerateVPRK(VPRK(tableau, args...; kwargs...))
DegenerateVPRK(method::Union{RKMethod,PRKMethod}, args...; kwargs...) = DegenerateVPRK(VPRK(method, args...; kwargs...))

tableau(method::DegenerateVPRK) = tableau(method.vprk)
nullvector(method::DegenerateVPRK) = nullvector(method.vprk)
hasnullvector(method::DegenerateVPRK) = hasnullvector(method.vprk)


function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::DegenerateVPRK; kwargs...)
    IntegratorVPRKdegenerate(problem, tableau(method), nullvector(method); kwargs...)
end

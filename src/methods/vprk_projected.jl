"""
Projected [`VPRK`](@ref) Method

A projected Variational Runge-Kutta Method consists of a VPRK method and a projection method:
```
ProjectedVPRK(vprk::VPRKMethod, projection::ProjectionMethod)
ProjectedVPRK(tableau::AbstractTableau, projection::ProjectionMethod, args...; kwargs...) = ProjectedVPRK(VPRK(tableau, args...; kwargs...), projection)
```
"""
struct ProjectedVPRK{VT,PT} <: VPRKMethod
    vprk::VT
    projection::PT

    function ProjectedVPRK(vprk::VT, projection::PT) where {VT <: VPRKMethod, PT <: ProjectionMethod}
        new{VT,PT}(vprk, projection)
    end
end

ProjectedVPRK(tableau::AbstractTableau, projection::ProjectionMethod, args...; kwargs...) = ProjectedVPRK(VPRK(tableau, args...; kwargs...), projection)

tableau(method::ProjectedVPRK) = tableau(method.vprk)
nullvector(method::ProjectedVPRK) = nullvector(method.vprk)
hasnullvector(method::ProjectedVPRK) = hasnullvector(method.vprk)


"[`VPRK`](@ref) integrator with projection on internal stages"
VPRKpInternal(method::VPRKMethod) = ProjectedVPRK(method, InternalStageProjection())

"[`VPRK`](@ref) integrator with Legendre projection"
VPRKpLegendre(method::VPRKMethod) = ProjectedVPRK(method, LegendreProjection())

"[`VPRK`](@ref) integrator with Midpoint projection"
VPRKpMidpoint(method::VPRKMethod) = ProjectedVPRK(method, MidpointProjection())

"[`VPRK`](@ref) integrator with projection on secondary constraint"
VPRKpSecondary(method::VPRKMethod) = ProjectedVPRK(method, SecondaryProjection())

"[`VPRK`](@ref) integrator with standard projection"
VPRKpStandard(method::VPRKMethod) = ProjectedVPRK(method, StandardProjection())

"[`VPRK`](@ref) integrator with symmetric projection"
VPRKpSymmetric(method::VPRKMethod) = ProjectedVPRK(method, SymmetricProjection())

"[`VPRK`](@ref) integrator with symplectic projection"
VPRKpSymplectic(method::VPRKMethod) = ProjectedVPRK(method, SymplecticProjection())

"[`VPRK`](@ref) integrator with variational projection"
VPRKpVariational(method::VPRKMethod) = ProjectedVPRK(method, VariationalProjection())

"[`VPRK`](@ref) integrator with variational projection on P"
VPRKpVariationalP(method::VPRKMethod) = ProjectedVPRK(method, VariationalProjectionOnP())

"[`VPRK`](@ref) integrator with variational projection on Q"
VPRKpVariationalQ(method::VPRKMethod) = ProjectedVPRK(method, VariationalProjectionOnQ())

VPRKpInternal(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), InternalStageProjection())
VPRKpLegendre(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), LegendreProjection())
VPRKpMidpoint(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), MidpointProjection())
VPRKpSecondary(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), SecondaryProjection())
VPRKpStandard(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), StandardProjection())
VPRKpSymmetric(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), SymmetricProjection())
VPRKpSymplectic(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), SymplecticProjection())
VPRKpVariational(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), VariationalProjection())
VPRKpVariationalP(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), VariationalProjectionOnP())
VPRKpVariationalQ(args...; kwargs...) = ProjectedVPRK(VPRK(args...; kwargs...), VariationalProjectionOnQ())


function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: InternalStageProjection}
    IntegratorVPRKpInternal(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: LegendreProjection}
    IntegratorVPRKpLegendre(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: MidpointProjection}
    IntegratorVPRKpMidpoint(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IDAEProblem,LDAEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SecondaryProjection}
    IntegratorVPRKpSecondary(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: StandardProjection}
    IntegratorVPRKpStandard(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SymmetricProjection}
    IntegratorVPRKpSymmetric(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SymplecticProjection}
    IntegratorVPRKpSymplectic(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjection}
    IntegratorVPRKpVariational(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjectionOnP}
    IntegratorVPRKpVariationalP(problem, tableau(method), nullvector(method); kwargs...)
end

function Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjectionOnQ}
    IntegratorVPRKpVariationalQ(problem, tableau(method), nullvector(method); kwargs...)
end

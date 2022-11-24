
struct ProjectedVPRK{VT,PT} <: VPRKMethod
    vprk::VT
    projection::PT

    function ProjectedVPRK(vprk::VT, projection::PT) where {VT <: VPRKMethod, PT <: ProjectionMethod}
        new{VT,PT}(vprk, projection)
    end
end

ProjectedVPRK(projection::ProjectionMethod, tableau::AbstractTableau, args...; kwargs...) = ProjectedVPRK(projection, VPRK(tableau, args...; kwargs...))

tableau(method::ProjectedVPRK) = tableau(method.vprk)
nullvector(method::ProjectedVPRK) = nullvector(method.vprk)
hasnullvector(method::ProjectedVPRK) = hasnullvector(method.vprk)

VPRKpInternal(method::VPRK) = ProjectedVPRK(method, InternalStageProjection())
VPRKpLegendre(method::VPRK) = ProjectedVPRK(method, LegendreProjection())
VPRKpMidpoint(method::VPRK) = ProjectedVPRK(method, MidpointProjection())
VPRKpSecondary(method::VPRK) = ProjectedVPRK(method, SecondaryProjection())
VPRKpStandard(method::VPRK) = ProjectedVPRK(method, StandardProjection())
VPRKpSymmetric(method::VPRK) = ProjectedVPRK(method, SymmetricProjection())
VPRKpSymplectic(method::VPRK) = ProjectedVPRK(method, SymplecticProjection())
VPRKpVariational(method::VPRK) = ProjectedVPRK(method, VariationalProjection())
VPRKpVariationalP(method::VPRK) = ProjectedVPRK(method, VariationalProjectionOnP())
VPRKpVariationalQ(method::VPRK) = ProjectedVPRK(method, VariationalProjectionOnQ())

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

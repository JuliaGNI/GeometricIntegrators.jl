
abstract type ProjectionMethod <: GeometricMethod end

struct NoProjection <: ProjectionMethod end

projection(::GeometricMethod) = NoProjection()


function _projection_weights(RU, RG, R∞ = 1)
    DT = promote_type(eltype(RU), eltype(RG), typeof(R∞))
    RU = Vector{DT}([RU[1], R∞ * RU[2]])
    RG = Vector{DT}([RG[1], R∞ * RG[2]])
    return (DT,RU,RG)
end

struct ProjectedMethod{PT <: ProjectionMethod, MT <: GeometricMethod} <: ProjectionMethod
    projection::PT
    method::MT

    function ProjectedMethod(proj::ProjectionMethod, method::GeometricMethod)
        _method = initmethod(method)
        new{typeof(proj), typeof(_method)}(proj, _method)
    end
end

projection(proj::ProjectedMethod) = proj.projection
Base.parent(proj::ProjectedMethod) = proj.method

initmethod(projection::ProjectedMethod) = projection
initmethod(projection::ProjectionMethod, method::GeometricMethod) = ProjectedMethod(projection, method)


struct InternalStageProjection <: ProjectionMethod end
struct LegendreProjection <: ProjectionMethod end
struct SecondaryProjection <: ProjectionMethod end
struct VariationalProjection <: ProjectionMethod end

InternalStageProjection(method::GeometricMethod) = ProjectedMethod(InternalStageProjection(), method)
LegendreProjection(method::GeometricMethod) = ProjectedMethod(LegendreProjection(), method)
SecondaryProjection(method::GeometricMethod) = ProjectedMethod(SecondaryProjection(), method)
VariationalProjection(method::GeometricMethod) = ProjectedMethod(VariationalProjection(), method)


struct MidpointProjection{DT} <: ProjectionMethod 
    RU::Vector{DT}
    RG::Vector{DT}

    function MidpointProjection(R∞ = 1)
        DT, RU, RG = _projection_weights([1//2,1//2], [1//2,1//2], R∞)
        new{DT}(RU, RG)
    end
end

MidpointProjection(method::GeometricMethod) = ProjectedMethod(MidpointProjection(), method)
MidpointProjection(method::Union{RKMethod,PRKMethod,VPRKMethod}) = ProjectedMethod(MidpointProjection(tableau(method).R∞), method)


struct SymmetricProjection{DT} <: ProjectionMethod 
    RU::Vector{DT}
    RG::Vector{DT}
    
    function SymmetricProjection(R∞ = 1)
        DT, RU, RG = _projection_weights([1//2,1//2], [1//2,1//2], R∞)
        new{DT}(RU, RG)
    end
end

SymmetricProjection(method::GeometricMethod) = ProjectedMethod(SymmetricProjection(), method)
SymmetricProjection(method::Union{RKMethod,PRKMethod,VPRKMethod}) = ProjectedMethod(SymmetricProjection(tableau(method).R∞), method)


struct StandardProjection{DT} <: ProjectionMethod 
    RU::Vector{DT}
    RG::Vector{DT}
    
    function StandardProjection(RU, RG, R∞ = 1)
        DT, RU, RG = _projection_weights(RU, RG, R∞)
        new{DT}(RU, RG)
    end
end

PostProjection(method::GeometricMethod) = ProjectedMethod(StandardProjection([0,1], [0,1]), method)
SymplecticProjection(method::Union{PRKMethod,VPRKMethod}) = ProjectedMethod(StandardProjection([0,1], [0,1], tableau(method).R∞), method)
VariationalProjectionOnP(method::GeometricMethod) = ProjectedMethod(StandardProjection([0,1], [1,0]), method)
VariationalProjectionOnQ(method::GeometricMethod) = ProjectedMethod(StandardProjection([1,0], [0,1]), method)

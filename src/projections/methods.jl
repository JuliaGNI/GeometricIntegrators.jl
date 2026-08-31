struct InternalStageProjection <: ProjectionMethod end
struct LegendreProjection <: ProjectionMethod end
struct SecondaryProjection <: ProjectionMethod end
struct VariationalProjection <: ProjectionMethod end

function InternalStageProjection(method::GeometricMethod)
    ProjectedMethod(InternalStageProjection(), method)
end
LegendreProjection(method::GeometricMethod) = ProjectedMethod(LegendreProjection(), method)
function SecondaryProjection(method::GeometricMethod)
    ProjectedMethod(SecondaryProjection(), method)
end
function VariationalProjection(method::GeometricMethod)
    ProjectedMethod(VariationalProjection(), method)
end

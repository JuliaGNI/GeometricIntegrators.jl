struct InternalStageProjection <: ProjectionMethod end
struct LegendreProjection <: ProjectionMethod end
struct SecondaryProjection <: ProjectionMethod end
struct VariationalProjection <: ProjectionMethod end

InternalStageProjection(method::GeometricMethod) = ProjectedMethod(InternalStageProjection(), method)
LegendreProjection(method::GeometricMethod) = ProjectedMethod(LegendreProjection(), method)
SecondaryProjection(method::GeometricMethod) = ProjectedMethod(SecondaryProjection(), method)
VariationalProjection(method::GeometricMethod) = ProjectedMethod(VariationalProjection(), method)

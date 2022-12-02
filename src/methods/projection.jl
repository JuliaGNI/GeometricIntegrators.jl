
abstract type ProjectionMethod <: Method end

struct InternalStageProjection <: ProjectionMethod end
struct LegendreProjection <: ProjectionMethod end
struct MidpointProjection <: ProjectionMethod end
struct SecondaryProjection <: ProjectionMethod end
struct StandardProjection <: ProjectionMethod end
struct SymmetricProjection <: ProjectionMethod end
struct SymplecticProjection <: ProjectionMethod end
struct VariationalProjection <: ProjectionMethod end
struct VariationalProjectionOnP <: ProjectionMethod end
struct VariationalProjectionOnQ <: ProjectionMethod end

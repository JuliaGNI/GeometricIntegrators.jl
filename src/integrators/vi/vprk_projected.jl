
"[`VPRK`](@ref) integrator with projection on internal stages"
VPRKpInternal(method::VPRKMethod) = InternalStageProjection(method)

"[`VPRK`](@ref) integrator with Legendre projection"
VPRKpLegendre(method::VPRKMethod) = LegendreProjection(method)

"[`VPRK`](@ref) integrator with Midpoint projection"
VPRKpMidpoint(method::VPRKMethod) = MidpointProjection(method)

"[`VPRK`](@ref) integrator with projection on secondary constraint"
VPRKpSecondary(method::VPRKMethod) = SecondaryProjection(method)

"[`VPRK`](@ref) integrator with standard projection"
VPRKpStandard(method::VPRKMethod) = PostProjection(method)

"[`VPRK`](@ref) integrator with symmetric projection"
VPRKpSymmetric(method::VPRKMethod) = SymmetricProjection(method)

"[`VPRK`](@ref) integrator with symplectic projection"
VPRKpSymplectic(method::VPRKMethod) = SymplecticProjection(method)

"[`VPRK`](@ref) integrator with variational projection"
VPRKpVariational(method::VPRKMethod) = VariationalProjection(method)

"[`VPRK`](@ref) integrator with variational projection on P"
VPRKpVariationalP(method::VPRKMethod) = VariationalProjectionOnP(method)

"[`VPRK`](@ref) integrator with variational projection on Q"
VPRKpVariationalQ(method::VPRKMethod) = VariationalProjectionOnQ(method)

VPRKpInternal(args...; kwargs...) = InternalStageProjection(VPRK(args...; kwargs...))
VPRKpLegendre(args...; kwargs...) = LegendreProjection(VPRK(args...; kwargs...))
VPRKpMidpoint(args...; kwargs...) = MidpointProjection(VPRK(args...; kwargs...))
VPRKpSecondary(args...; kwargs...) = SecondaryProjection(VPRK(args...; kwargs...))
VPRKpStandard(args...; kwargs...) = PostProjection(VPRK(args...; kwargs...))
VPRKpSymmetric(args...; kwargs...) = SymmetricProjection(VPRK(args...; kwargs...))
VPRKpSymplectic(args...; kwargs...) = SymplecticProjection(VPRK(args...; kwargs...))
VPRKpVariational(args...; kwargs...) = VariationalProjection(VPRK(args...; kwargs...))
VPRKpVariationalP(args...; kwargs...) = VariationalProjectionOnP(VPRK(args...; kwargs...))
VPRKpVariationalQ(args...; kwargs...) = VariationalProjectionOnQ(VPRK(args...; kwargs...))

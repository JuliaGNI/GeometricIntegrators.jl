
# function check_symplecticity end
function symplecticity_conditions end



# GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::RKMethod; kwargs...) = GeometricIntegrator(problem, Methods.tableau(method); kwargs...)
# GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::PRKMethod; kwargs...) = GeometricIntegrator(problem, Methods.tableau(method); kwargs...)

# GeometricIntegrator(problem::LODEProblem, method::FLRK; kwargs...) = IntegratorFLRK(problem, Methods.tableau(method); kwargs...)

# GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::VPRKMethod; kwargs...) = IntegratorVPRK(problem, Methods.tableau(method), nullvector(method); kwargs...)

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: InternalStageProjection}
#     IntegratorVPRKpInternal(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: LegendreProjection}
#     IntegratorVPRKpLegendre(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: MidpointProjection}
#     IntegratorVPRKpMidpoint(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IDAEProblem,LDAEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SecondaryProjection}
#     IntegratorVPRKpSecondary(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: StandardProjection}
#     IntegratorVPRKpStandard(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SymmetricProjection}
#     IntegratorVPRKpSymmetric(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SymplecticProjection}
#     IntegratorVPRKpSymplectic(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjection}
#     IntegratorVPRKpVariational(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjectionOnP}
#     IntegratorVPRKpVariationalP(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjectionOnQ}
#     IntegratorVPRKpVariationalQ(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function GeometricIntegrator(problem::Union{IODEProblem,LODEProblem}, method::DegenerateVPRK; kwargs...)
#     IntegratorVPRKdegenerate(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# GeometricIntegrator(problem::LODEProblem, method::DVIA; kwargs...) = IntegratorDVIA(problem; kwargs...)
# GeometricIntegrator(problem::LODEProblem, method::DVIB; kwargs...) = IntegratorDVIB(problem; kwargs...)
# GeometricIntegrator(problem::LODEProblem, method::CMDVI; kwargs...) = IntegratorCMDVI(problem; kwargs...)
# GeometricIntegrator(problem::LODEProblem, method::CTDVI; kwargs...) = IntegratorCTDVI(problem; kwargs...)

# GeometricIntegrator(problem::SODEProblem, method::SplittingMethod; kwargs...) = GeometricIntegrator(problem, Methods.tableau(method); kwargs...)
# GeometricIntegratorComposition(problem::SODEProblem, method::SplittingMethod; kwargs...) = IntegratorComposition(problem, Methods.tableau(method); kwargs...)
# GeometricIntegratorComposition(problem::SODEProblem, integrators::Tuple, method::SplittingMethod; kwargs...) = IntegratorComposition(problem, integrators, Methods.tableau(method); kwargs...)

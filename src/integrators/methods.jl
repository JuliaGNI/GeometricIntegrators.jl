
# function check_symplecticity end
function symplecticity_conditions end



# AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::RKMethod; kwargs...) = AbstractIntegrator(problem, Methods.tableau(method); kwargs...)
# AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::PRKMethod; kwargs...) = AbstractIntegrator(problem, Methods.tableau(method); kwargs...)

AbstractIntegrator(problem::LODEProblem, method::FLRK; kwargs...) = IntegratorFLRK(problem, Methods.tableau(method); kwargs...)

# AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::VPRKMethod; kwargs...) = IntegratorVPRK(problem, Methods.tableau(method), nullvector(method); kwargs...)

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: InternalStageProjection}
#     IntegratorVPRKpInternal(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: LegendreProjection}
#     IntegratorVPRKpLegendre(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: MidpointProjection}
#     IntegratorVPRKpMidpoint(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IDAEProblem,LDAEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SecondaryProjection}
#     IntegratorVPRKpSecondary(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: StandardProjection}
#     IntegratorVPRKpStandard(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SymmetricProjection}
#     IntegratorVPRKpSymmetric(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: SymplecticProjection}
#     IntegratorVPRKpSymplectic(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjection}
#     IntegratorVPRKpVariational(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjectionOnP}
#     IntegratorVPRKpVariationalP(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::ProjectedVPRK{VT,PT}; kwargs...) where {VT, PT <: VariationalProjectionOnQ}
#     IntegratorVPRKpVariationalQ(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, method::DegenerateVPRK; kwargs...)
#     IntegratorVPRKdegenerate(problem, Methods.tableau(method), nullvector(method); kwargs...)
# end

# AbstractIntegrator(problem::LODEProblem, method::DVIA; kwargs...) = IntegratorDVIA(problem; kwargs...)
# AbstractIntegrator(problem::LODEProblem, method::DVIB; kwargs...) = IntegratorDVIB(problem; kwargs...)
# AbstractIntegrator(problem::LODEProblem, method::CMDVI; kwargs...) = IntegratorCMDVI(problem; kwargs...)
# AbstractIntegrator(problem::LODEProblem, method::CTDVI; kwargs...) = IntegratorCTDVI(problem; kwargs...)

# AbstractIntegrator(problem::SODEProblem, method::SplittingMethod; kwargs...) = AbstractIntegrator(problem, Methods.tableau(method); kwargs...)
# AbstractIntegratorComposition(problem::SODEProblem, method::SplittingMethod; kwargs...) = IntegratorComposition(problem, Methods.tableau(method); kwargs...)
# AbstractIntegratorComposition(problem::SODEProblem, integrators::Tuple, method::SplittingMethod; kwargs...) = IntegratorComposition(problem, integrators, Methods.tableau(method); kwargs...)

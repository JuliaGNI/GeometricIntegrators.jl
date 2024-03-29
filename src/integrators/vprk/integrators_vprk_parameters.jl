
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct AbstractParametersVPRK{IT, DT, TT, D, S, ET, PT, NV} <: Parameters{DT,TT}
    equ::ET
    tab::PartitionedTableau{TT}
    nullvec::NV
    Δt::TT

    pparams::PT

    t̄::TT
    q̄::Vector{DT}
    p̄::Vector{DT}
    v̄::Vector{DT}

    function AbstractParametersVPRK{IT,DT,D}(equs::ET, tab::PartitionedTableau{TT}, nullvec::NV, Δt::TT, pparams::PT=NamedTuple()) where {IT, DT, TT, D, ET <: NamedTuple, PT <: NamedTuple, NV <: Union{AbstractArray,Nothing}}
        new{IT, DT, TT, D, tab.s, ET, PT, NV}(equs, tab, nullvec, Δt, pparams, zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function update_params!(params::AbstractParametersVPRK, sol::Union{SolutionStepPODE, SolutionStepPDAE})
    # set time for nonlinear solver and copy previous solution
    solstep.t̄  = sol.t
    solstep.q̄ .= sol.q
    solstep.p̄ .= sol.p
    params.v̄ .= sol.v
end

@inline equation(int::GeometricIntegratorVPRK, i::Symbol) = parameters(int).equ[i]
@inline equations(int::GeometricIntegratorVPRK) = parameters(int).equ
@inline GeometricBase.timestep(int::GeometricIntegratorVPRK) = parameters(int).Δt
@inline GeometricBase.tableau(int::GeometricIntegratorVPRK) = parameters(int).tab


function Integrators.IntegratorCache{ST}(params::AbstractParametersVPRK{IT,DT,TT,D,S}; kwargs...) where {IT,ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(true; kwargs...)
end

@inline Integrators.CacheType(ST, params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,DT,TT,D,S} = IntegratorCacheVPRK{ST,D,S}


"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct AbstractParametersVPRK{IT, DT, TT, D, S, ET <: NamedTuple, PT <: NamedTuple} <: Parameters{DT,TT}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    pparams::PT

    t̅::TT
    q̅::Vector{DT}
    p̅::Vector{DT}

    function AbstractParametersVPRK{IT,DT,D}(equs::ET, tab::TableauVPRK{TT}, Δt::TT, pparams::PT=NamedTuple()) where {IT, DT, TT, D, ET <: NamedTuple, PT <: NamedTuple}
        new{IT, DT, TT, D, tab.s, ET, PT}(equs, tab, Δt, pparams, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end

function update_params!(params::AbstractParametersVPRK, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t̅  = sol.t
    params.q̅ .= sol.q
    params.p̅ .= sol.p
end

@inline equation(int::AbstractIntegratorVPRK, i::Symbol) = parameters(int).equ[i]
@inline equations(int::AbstractIntegratorVPRK) = parameters(int).equ
@inline timestep(int::AbstractIntegratorVPRK) = parameters(int).Δt
@inline tableau(int::AbstractIntegratorVPRK) = parameters(int).tab


function Integrators.IntegratorCache(params::AbstractParametersVPRK{IT,DT,TT,D,S}; kwargs...) where {IT,DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(true; kwargs...)
end

function Integrators.IntegratorCache{ST}(params::AbstractParametersVPRK{IT,DT,TT,D,S}; kwargs...) where {IT,ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(true; kwargs...)
end

@inline Integrators.CacheType(ST, params::AbstractParametersVPRK{IT,DT,TT,D,S}) where {IT,DT,TT,D,S} = IntegratorCacheVPRK{ST,D,S}

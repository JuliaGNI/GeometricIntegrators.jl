
"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods."
mutable struct AbstractParametersSPARK{IT, DT, TT, D, S, R, P, equType <: NamedTuple, tabType <: AbstractTableau} <: Parameters{DT,TT}
    equs::equType
    tab::tabType
    Δt::TT

    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    function AbstractParametersSPARK{IT,DT,D,S,R,P}(equs::equType, tab::tabType, Δt::TT) where {IT,DT,TT,D,S,R,P,equType,tabType}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{IT, DT, TT, D, S, R, P, equType, tabType}(equs, tab, Δt, zero(TT), q, p, λ)
    end

    function AbstractParametersSPARK{IT,DT,D}(equs::NamedTuple, tab::AbstractTableauSPARK, Δt::Number) where {IT,DT,D}
        AbstractParametersSPARK{IT, DT, D, tab.s, tab.r, tab.ρ}(equs, tab, Δt)
    end

    function AbstractParametersSPARK{IT,DT,D}(equs::NamedTuple, tab::Union{AbstractTableau,Tableau}, Δt::Number) where {IT,DT,D}
        AbstractParametersSPARK{IT, DT, D, tab.s, tab.r, 0}(equs, tab, Δt)
    end
end

function update_params!(params::AbstractParametersSPARK, sol::AtomicSolutionPDAE)
    # set time for nonlinear solver and copy previous solution
    params.t  = sol.t
    params.q .= sol.q
    params.p .= sol.p
    params.λ .= sol.λ
end

@inline equation(int::AbstractIntegratorSPARK, i::Symbol) = parameters(int).equs[i]
@inline equations(int::AbstractIntegratorSPARK) = parameters(int).equs
@inline timestep(int::AbstractIntegratorSPARK) = parameters(int).Δt
@inline tableau(int::AbstractIntegratorSPARK) = parameters(int).tab


function Integrators.IntegratorCache{ST}(params::AbstractParametersSPARK{IT,DT,TT,D,S,R}; kwargs...) where {IT,ST,DT,TT,D,S,R}
    IntegratorCacheSPARK{ST,D,S,R}(; kwargs...)
end

@inline Integrators.CacheType(ST, params::AbstractParametersSPARK{IT,DT,TT,D,S,R}) where {IT,DT,TT,D,S,R} = IntegratorCacheSPARK{ST,D,S,R}

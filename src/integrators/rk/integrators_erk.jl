
"Parameters for right-hand side function of explicit Runge-Kutta methods."
struct ParametersERK{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::Tableau{TT}
    Δt::TT

    function ParametersERK{DT,D}(equs::ET, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt)
    end
end


"Explicit Runge-Kutta integrator cache."
struct IntegratorCacheERK{DT,D,S} <: ODEIntegratorCache{DT,D}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}

    function IntegratorCacheERK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        new(Q, V)
    end
end

function IntegratorCache{ST}(params::ParametersERK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheERK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersERK{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheERK{ST,D,S}


"Explicit Runge-Kutta integrator."
struct IntegratorERK{DT, TT, D, S, ET} <: AbstractIntegratorRK{DT,TT}
    params::ParametersERK{DT,TT,D,S,ET}
    caches::CacheDict{ParametersERK{DT,TT,D,S,ET}}


    function IntegratorERK(params::ParametersERK{DT,TT,D,S,ET}, caches) where {DT,TT,D,S,ET}
        new{DT, TT, D, S, ET}(params, caches)
    end

    function IntegratorERK{DT,D}(equations::ET, tableau::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        # get number of stages
        S = tableau.s

        # check if tableau is explicit
        @assert isexplicit(tableau)

        # create params
        params = ParametersERK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create integrator
        IntegratorERK(params, caches)
    end

    function IntegratorERK{DT,D}(v::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorERK{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    end

    function IntegratorERK{DT,D}(v::Function, h::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorERK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorERK(equation::ODE{DT}, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorERK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorERK{DT,TT,D,S}) where {DT,TT,D,S} = D


function integrate_step!(int::IntegratorERK{DT,TT}, sol::AtomicSolutionODE{DT,TT},
                         cache::IntegratorCacheERK{DT}=int.caches[DT]) where {DT,TT}
    # temporary variables
    local tᵢ::TT
    local yᵢ::DT

    # reset atomic solution
    reset!(sol, timestep(int))

    # compute internal stages
    for i in eachstage(int)
        for k in eachindex(cache.Q[i], cache.V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += tableau(int).a[i,j] * cache.V[j][k]
            end
            cache.Q[i][k] = sol.q̄[k] + timestep(int) * yᵢ
        end
        tᵢ = sol.t̄ + timestep(int) * tableau(int).c[i]
        equations(int)[:v](tᵢ, cache.Q[i], cache.V[i])
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).b, tableau(int).b̂, timestep(int))
end

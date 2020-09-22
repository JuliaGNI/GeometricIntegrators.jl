
"Holds the tableau of an explicit Runge-Kutta method."
struct TableauERK{T} <: AbstractTableauERK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauERK{T}(q) where {T}
        @assert q.c[1] == 0
        @assert istrilstrict(q.a)
        @assert !(q.s==1 && q.a[1,1] ≠ 0)
        new(q.name, q.o, q.s, q)
    end
end

function TableauERK(q::CoefficientsRK{T}) where {T}
    TableauERK{T}(q)
end

function TableauERK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauERK{T}(CoefficientsRK(name, order, a, b, c))
end

"Read explicit Runge-Kutta tableau from file."
function readTableauERKFromFile(dir::AbstractString, name::AbstractString)
    file = string(dir, "/", name, ".tsv")

    o, s, T = readTableauRKHeaderFromFile(file)

    # TODO Read data in original format (e.g., Rational).
    #      For this we need to save tableaus as jld or hdf5.
#    tab_array = readdlm(file, T)
    tab_array = readdlm(file, comments=true)

    if s == 0
        s = size(tab_array, 1)
    end

    @assert s == size(tab_array, 1) == size(tab_array, 2)-2

    c = tab_array[1:s, 1]
    b = tab_array[1:s, 2]
    a = tab_array[1:s, 3:s+2]

    get_config(:verbosity) > 1 ? @info("Reading explicit Runge-Kutta tableau $(name) with $(s) stages and order $(o) from file\n$(file)") : nothing

    TableauERK(Symbol(name), o, a, b, c)
end


"Parameters for right-hand side function of explicit Runge-Kutta methods."
struct ParametersERK{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::TableauERK{TT}
    Δt::TT

    function ParametersERK{DT,D}(equs::ET, tab::TableauERK{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
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
struct IntegratorERK{DT, TT, D, S, ET} <: IntegratorRK{DT,TT}
    params::ParametersERK{DT,TT,D,S,ET}
    caches::CacheDict{ParametersERK{DT,TT,D,S,ET}}


    function IntegratorERK(params::ParametersERK{DT,TT,D,S,ET}, caches) where {DT,TT,D,S,ET}
        new{DT, TT, D, S, ET}(params, caches)
    end

    function IntegratorERK{DT,D}(equations::ET, tableau::TableauERK{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersERK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create integrator
        IntegratorERK(params, caches)
    end

    function IntegratorERK{DT,D}(v::Function, tableau::TableauERK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorERK{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    end

    function IntegratorERK{DT,D}(v::Function, h::Function, tableau::TableauERK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorERK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorERK(equation::ODE{DT,TT}, tableau::TableauERK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorERK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorERK{DT,TT,D,S}) where {DT,TT,D,S} = D


"Integrate ODE with explicit Runge-Kutta integrator."
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
                yᵢ += tableau(int).q.a[i,j] * cache.V[j][k]
            end
            cache.Q[i][k] = sol.q̅[k] + timestep(int) * yᵢ
        end
        tᵢ = sol.t̅ + timestep(int) * tableau(int).q.c[i]
        equations(int)[:v](tᵢ, cache.Q[i], cache.V[i])
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
end

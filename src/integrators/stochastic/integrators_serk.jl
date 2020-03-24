@doc raw"""
Holds the tableau of a stochastic explicit Runge-Kutta method.

Order of the tableau is not included, because unlike in the deterministic
setting, it depends on the properties of the noise (e.g., the dimension of
the Wiener process and the commutativity properties of the diffusion matrix).

Orders stored in `qdrift`, `qdiff` and `qdiff2` are understood as the classical orders of these methods.
"""
struct TableauSERK{T} <: AbstractTableauERK{T}
    name::Symbol
    s::Int

    qdrift::CoefficientsRK{T}
    qdiff::CoefficientsRK{T}
    qdiff2::CoefficientsRK{T}

    function TableauSERK{T}(name, qdrift, qdiff, qdiff2) where {T}
        @assert qdrift.s == qdiff.s == qdiff2.s
        @assert qdrift.c[1] == qdiff.c[1] == qdiff2.c[1] == 0.
        @assert istrilstrict(qdrift.a)
        @assert istrilstrict(qdiff.a)
        @assert istrilstrict(qdiff2.a)
        @assert !(qdrift.s==1. && qdrift.a[1,1] ≠ 0.)
        @assert !(qdiff.s==1. && qdiff.a[1,1] ≠ 0.)
        @assert !(qdiff2.s==1. && qdiff2.a[1,1] ≠ 0.)
        new(name, qdrift.s, qdrift, qdiff, qdiff2)
    end
end

function TableauSERK(name, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}, qdiff2::CoefficientsRK{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, qdiff2)
end

function TableauSERK(name, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, CoefficientsRK(T, :NULL, 0, zero(qdrift.a), zero(qdrift.b), zero(qdrift.c)) )
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T},
                                    order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T},
                                    order_diff2::Int, a_diff2::Matrix{T}, b_diff2::Vector{T}, c_diff2::Vector{T}) where {T}
    TableauSERK{T}(name, CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff),
                    CoefficientsRK(name, order_diff2, a_diff2, b_diff2, c_diff2))
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T},
                                    order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T}) where {T}
    TableauSERK{T}(name, CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff),
                    CoefficientsRK(:NULL, 0, zero(a_drift), zero(b_drift), zero(c_drift)))
end


"Parameters for stochastic explicit Runge-Kutta methods."
struct ParametersSERK{DT, TT, D, M, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::TableauSERK{TT}
    Δt::TT

    function ParametersSERK{DT,D,M}(equs::ET, tab::TableauSERK{TT}, Δt::TT) where {DT, TT, D, M, ET <: NamedTuple}
        new{DT, TT, D, M, tab.s, ET}(equs, tab, Δt)
    end
end


"Stochastic Explicit Runge-Kutta integrator cache."
struct IntegratorCacheSERK{DT,D,M,S} <: SDEIntegratorCache{DT,D,M}
    Δy::Vector{DT}

    Q::Vector{Vector{DT}}     # Q[j][k] - the k-th component of the j-th internal stage
    V::Vector{Vector{DT}}     # V[j][k] - the k-th component of v(Q[j])
    B::Vector{Matrix{DT}}     # B[j]    - the diffusion matrix B(Q[j])

    function IntegratorCacheSERK{DT,D,M,S}() where {DT,D,M,S}
        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        B = create_internal_stage_matrix(DT, D, M, S)

        new(zeros(DT,M), Q, V, B)
    end
end

function Integrators.IntegratorCache{ST}(params::ParametersSERK{DT,TT,D,M,S}; kwargs...) where {ST,DT,TT,D,M,S}
    IntegratorCacheSERK{ST,D,M,S}(; kwargs...)
end

@inline Integrators.CacheType(ST, params::ParametersSERK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = IntegratorCacheSERK{ST,D,M,S}


"Stochastic Explicit Runge-Kutta integrator."
struct IntegratorSERK{DT, TT, D, M, S, ET <: NamedTuple} <: StochasticIntegratorRK{DT,TT,D,M,S}
    params::ParametersSERK{DT,TT,D,M,S,ET}
    caches::CacheDict{ParametersSERK{DT,TT,D,M,S,ET}}

    function IntegratorSERK(params::ParametersSERK{DT,TT,D,M,S,ET}, caches) where {DT,TT,D,M,S,ET}
        new{DT, TT, D, M, S, ET}(params, caches)
    end

    function IntegratorSERK{DT,D,M}(equations::ET, tableau::TableauSERK{TT}, Δt::TT) where {DT, TT, D, M, ET <: NamedTuple}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersSERK{DT,D,M}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create integrator
        IntegratorSERK(params, caches)
    end

    function IntegratorSERK(equation::SDE{DT,TT}, tableau::TableauSERK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorSERK{DT, ndims(equation), equation.m}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


"Integrate SDE with explicit Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorSERK{DT,TT}, sol::AtomicSolutionSDE{DT,TT},
                                     cache::IntegratorCacheSERK{DT}=int.caches[DT]) where {DT,TT}
    local tᵢ::TT
    local ydrift::DT

    # reset cache
    reset!(sol, timestep(int))

    for i in eachstage(int)
        for k in eachindex(cache.Q[i])
            # contribution from the drift part
            ydrift = 0
            for j = 1:i-1
                ydrift += tableau(int).qdrift.a[i,j] * cache.V[j][k]
            end

            # ΔW contribution from the diffusion part
            cache.Δy .= 0
            for j = 1:i-1
                for l in eachnoise(int)
                    cache.Δy[l] += tableau(int).qdiff.a[i,j] * cache.B[j][k,l]
                end
            end

            cache.Q[i][k] = sol.q̅[k] + timestep(int) * ydrift + dot(cache.Δy, sol.ΔW)

            # ΔZ contribution from the diffusion part
            if tableau(int).qdiff2.name ≠ :NULL
                cache.Δy .= 0
                for j = 1:i-1
                    for l in eachnoise(int)
                        cache.Δy[l] += tableau(int).qdiff2.a[i,j] * cache.B[j][k,l]
                    end
                end

                cache.Q[i][k] += dot(cache.Δy, sol.ΔZ)/timestep(int)
            end

        end

        tᵢ = sol.t̅ + timestep(int) * tableau(int).qdrift.c[i]
        equation(int, :v)(tᵢ, cache.Q[i], cache.V[i])
        equation(int, :B)(tᵢ, cache.Q[i], cache.B[i])
    end

    # compute final update
    if tableau(int).qdiff2.name == :NULL
        update_solution!(sol, cache.V, cache.B, tableau(int).qdrift.b, tableau(int).qdiff.b, timestep(int), sol.ΔW, cache.Δy)
        # update_solution!(sol, cache.V, cache.B, tableau(int).qdrift.b̂, tableau(int).qdiff.b̂, timestep(int), int.ΔW, cache.Δy)
    else
        update_solution!(sol, cache.V, cache.B, tableau(int).qdrift.b, tableau(int).qdiff.b, tableau(int).qdiff2.b, timestep(int), sol.ΔW, sol.ΔZ, cache.Δy)
        # update_solution!(sol, cache.V, cache.B, tableau(int).qdrift.b̂, tableau(int).qdiff.b̂, tableau(int).qdiff2.b̂, timestep(int), sol.ΔW, sol.ΔZ, cache.Δy)
    end
end

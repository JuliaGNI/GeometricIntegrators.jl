
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

#    run(`cat $file`)

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

    @info "Reading explicit Runge-Kutta tableau $(name) with $(s) stages and order $(o) from file\n$(file)"

    TableauERK(Symbol(name), o, a, b, c)
end



"Explicit Runge-Kutta integrator."
struct IntegratorERK{DT,TT,FT} <: DeterministicIntegrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauERK{TT}
    Δt::TT

    function IntegratorERK{DT,TT,FT}(equation, tableau, Δt) where {DT,TT,FT}
        new(equation, tableau, Δt)
    end
end

function IntegratorERK(equation::ODE{DT,TT,FT}, tableau::TableauERK{TT}, Δt::TT) where {DT,TT,FT}
    IntegratorERK{DT,TT,FT}(equation, tableau, Δt)
end

equation(int::IntegratorERK) = int.equation
timestep(int::IntegratorERK) = int.Δt


"Explicit Runge-Kutta integrator cache."
mutable struct IntegratorCacheERK{DT,TT,D,S} <: ODEIntegratorCache{DT,D,S}
    t::TT
    t̅::TT
    q::Vector{DT}
    q̅::Vector{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}

    function IntegratorCacheERK{DT,TT,D,S}() where {DT,TT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        new(zero(TT), zero(TT), zeros(DT,D), zeros(DT,D), Q, V)
    end
end

function create_integrator_cache(int::IntegratorERK{DT,TT}) where {DT,TT}
    IntegratorCacheERK{DT, TT, int.equation.d, int.tableau.s}()
end

function CommonFunctions.reset!(cache::IntegratorCacheERK{DT,TT}, Δt::TT) where {DT,TT}
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.t += Δt
end

function CommonFunctions.get_solution(cache::IntegratorCacheERK)
    (cache.t, cache.q)
end

function CommonFunctions.set_solution!(cache::IntegratorCacheERK, sol)
    t, q = sol
    cache.t  = t
    cache.q .= q
end


"Integrate ODE with explicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorERK{DT,TT}, cache::IntegratorCacheERK{DT,TT}) where {DT,TT}
    # @assert m ≥ 1
    # @assert m ≤ sol.ni
    #
    # @assert n ≥ 1
    # @assert n ≤ sol.ntime
    #
    # t = sol.t[0] + (n-1)*int.Δt

    local tᵢ::TT
    local yᵢ::DT

    # reset cache
    reset!(cache, int.Δt)

    # compute internal stages
    for i in 1:int.tableau.q.s
        @inbounds for k in eachindex(cache.Q[i], cache.V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += int.tableau.q.a[i,j] * cache.V[j][k]
            end
            cache.Q[i][k] = cache.q̅[k] + int.Δt * yᵢ
        end
        tᵢ = cache.t̅ + int.Δt * int.tableau.q.c[i]
        int.equation.v(tᵢ, cache.Q[i], cache.V[i])
    end

    # compute final update
    update_solution!(cache.q, cache.V, int.tableau.q.b, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(cache.q, int.equation.periodicity)

    # copy to solution
    # copy_solution!(sol, int.q, n, m)
end

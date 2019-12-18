
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
struct ParametersERK{DT, TT, ET <: ODE{DT,TT}, D, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauERK{TT}
    Δt::TT
end

function ParametersERK(equ::ET, tab::TableauERK{TT}, Δt::TT) where {DT, TT, ET <: ODE{DT,TT}}
    ParametersERK{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt)
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


"Explicit Runge-Kutta integrator."
struct IntegratorERK{DT, TT, PT <: ParametersERK{DT,TT}, D, S} <: IntegratorRK{DT,TT}
    params::PT
    cache::IntegratorCacheERK{DT,D,S}

    function IntegratorERK(equation::ODE{DT,TT,FT}, tableau::TableauERK{TT}, Δt::TT) where {DT,TT,FT}
        D = equation.d
        S = tableau.s

        # create params
        params = ParametersERK(equation, tableau, Δt)

        # create cache
        cache = IntegratorCacheERK{DT,D,S}()

        # create integrators
        new{DT, TT, typeof(params), D, S}(params, cache)
    end
end


"Integrate ODE with explicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorERK{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    local tᵢ::TT
    local yᵢ::DT

    # reset cache
    reset!(sol, timestep(int))

    # compute internal stages
    for i in eachstage(int)
        @inbounds for k in eachindex(int.cache.Q[i], int.cache.V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += tableau(int).q.a[i,j] * int.cache.V[j][k]
            end
            int.cache.Q[i][k] = sol.q̅[k] + timestep(int) * yᵢ
        end
        tᵢ = sol.t̅ + timestep(int) * tableau(int).q.c[i]
        equation(int).v(tᵢ, int.cache.Q[i], int.cache.V[i])
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
end

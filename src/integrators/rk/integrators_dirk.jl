
"Holds the tableau of a diagonally implicit Runge-Kutta method."
struct TableauDIRK{T} <: AbstractTableauIRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauDIRK{T}(q) where {T}
        @assert istril(q.a)
        @assert !(q.s==1 && q.a[1,1] ≠ 0)

        if q.s > 1 && istrilstrict(q.a)
            @warn "Initializing TableauDIRK with explicit tableau $(q.name).\nYou might want to use TableauERK instead."
        end

        new(q.name, q.o, q.s, q)
    end
end

function TableauDIRK(q::CoefficientsRK{T}) where {T}
    TableauDIRK{T}(q)
end

function TableauDIRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauDIRK{T}(CoefficientsRK(name, order, a, b, c))
end

# TODO function readTableauDIRKFromFile(dir::AbstractString, name::AbstractString)


"Parameters for right-hand side function of diagonally implicit Runge-Kutta methods."
mutable struct ParametersDIRK{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::TableauDIRK{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    V::Vector{Vector{DT}}

    function ParametersDIRK{DT,D}(equs::ET, tab::TableauDIRK{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt, zero(TT), zeros(DT,D), create_internal_stage_vector(DT, D, tab.s))
    end
end


"""
Diagonally implicit Runge-Kutta integrator cache.
"""
struct IntegratorCacheDIRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}
    y::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCacheDIRK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Y)
    end
end


"Diagonally implicit Runge-Kutta integrator."
struct IntegratorDIRK{DT, TT, D, S, PT <: ParametersDIRK{DT,TT},
                                    ST,# <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{DT,TT}} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheDIRK{DT,D,S}

    function IntegratorDIRK(params::ParametersDIRK{DT,TT,D,S}, solver::ST, iguess::IT, cache) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, cache)
    end

    function IntegratorDIRK{DT,D}(equations::NamedTuple, tableau::TableauDIRK{TT}, Δt::TT) where {DT, TT, D}
        # get number of stages
        S = tableau.q.s

        # create params
        params = ParametersDIRK{DT,D}(equations, tableau, Δt)

        # create solvers
        function create_nonlinear_solver(DT, N, params, i)
            # create solution vector for nonlinear solver
            x = zeros(DT, N)

            # create wrapper function f(x,b) that calls `function_stages!(x, b, params)`
            # with the appropriate `params`
            f = (x,b) -> function_stages!(x, b, params, i)

            # create nonlinear solver with solver type obtained from config dictionary
            s = get_config(:nls_solver)(x, f)
        end

        solvers = [create_nonlinear_solver(DT, D, params, i) for i in 1:S]

        # create initial guess
        iguess = InitialGuessODE{DT,D}(get_config(:ig_interpolation), equations[:v], Δt)

        # create cache
        cache = IntegratorCacheDIRK{DT,D,S}()

        # create integrator
        IntegratorDIRK(params, solvers, iguess, cache)
    end

    function IntegratorDIRK{DT,D}(v::Function, tableau::TableauDIRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorDIRK{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    end

    function IntegratorDIRK{DT,D}(v::Function, h::Function, tableau::TableauDIRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorDIRK{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    end

    function IntegratorDIRK(equation::ODE{DT,TT}, tableau::TableauDIRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorDIRK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(int::IntegratorDIRK{DT,TT,D,S}) where {DT,TT,D,S} = D
@inline has_initial_guess(int::IntegratorDIRK) = true


"Initialise initial guess"
function initialize!(int::IntegratorDIRK, cache::IntegratorCacheDIRK)
    # initialise initial guess
    cache.t̅ = cache.t - timestep(int)

    int.params.equs[:v](cache.t, cache.q, cache.v)

    initialize!(int.iguess, cache.t, cache.q, cache.v, cache.t̅, cache.q̅, cache.v̅)
end


"Compute initial guess for internal stages."
function initial_guess!(int::IntegratorDIRK, sol::AtomicSolutionODE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.v, sol.q̅, sol.v̅, int.cache.q̃, int.cache.ṽ, int.params.tab.q.c[i])
        for k in eachindex(int.cache.V[i], int.cache.ṽ)
            int.cache.V[i][k] = int.cache.ṽ[k]
        end
    end
    for i in eachstage(int)
        for k in eachdim(int)
            int.solver[i].x[k] = 0
            for j in eachstage(int)
                int.solver[i].x[k] += int.params.tab.q.a[i,j] * int.cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{ST}, V::Vector{ST}, Y::Vector{ST},
                                    params::ParametersDIRK{DT,TT,D,S}, i::Int) where {ST,DT,TT,D,S}

    local tᵢ::TT

    @assert D == length(Q) == length(V) == length(Y)

    # copy x to Y and compute Q = q + Δt Y
    for k in 1:D
        Y[k] = x[k]
        Q[k] = params.q[k] + params.Δt * Y[k]
    end

    # compute V = v(Q)
    tᵢ = params.t + params.Δt * params.tab.q.c[i]
    params.equs[:v](tᵢ, Q, V)
end


"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDIRK{DT,TT,D,S}, i::Int) where {ST,DT,TT,D,S}
    local Q = zeros(ST,D)
    local V = zeros(ST,D)
    local Y = zeros(ST,D)

    local y1::ST
    local y2::ST

    compute_stages!(x, Q, V, Y, params, i)

    # compute b = - (Y-AV)
    for k in 1:D
        y1 = params.tab.q.a[i,i] * V[k]
        y2 = params.tab.q.â[i,i] * V[k]
        for j in 1:i-1
            y1 += params.tab.q.a[i,j] * params.V[j][k]
            y2 += params.tab.q.â[i,j] * params.V[j][k]
        end
        b[k] = - Y[k] + (y1 + y2)
    end
end


"Integrate ODE with diagonally implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorDIRK{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT,N}
    int.params.t  = sol.t
    int.params.q .= sol.q

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    # consecutively solve for all stages
    for i in eachstage(int)
        # call nonlinear solver
        solve!(int.solver[i])

        # print solver status
        print_solver_status(int.solver[i].status, int.solver[i].params)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.solver[i].status, int.solver[i].params)

        # compute vector field at internal stages
        compute_stages!(int.solver[i].x, int.cache.Q[i], int.cache.V[i], int.cache.Y[i], int.params, i)

        # copy velocity field to parameters
        int.params.V[i] .= int.cache.V[i]
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)

    # update vector field for initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end

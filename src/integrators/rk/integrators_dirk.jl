
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
mutable struct ParametersDIRK{DT, TT, ET <: ODE{DT,TT}, D, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauDIRK{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
end

function ParametersDIRK(equ::ET, tab::TableauDIRK{TT}, Δt::TT) where {DT, TT, ET <: ODE{DT,TT}}
    ParametersDIRK{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, 0, zeros(DT, equ.d))
end


struct NonlinearFunctionCacheDIRK{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    v::Vector{DT}
    y::Vector{DT}

    function NonlinearFunctionCacheDIRK{DT}(D,S) where {DT}
        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)

        # create velocity and update vector
        v = zeros(DT,D)
        y = zeros(DT,D)

        new(Q, V, Y, v, y)
    end
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                                    params::ParametersDIRK{DT,TT,ET,D,S}, i::Int) where {ST,DT,TT,ET,D,S}

    local tᵢ::TT

    @assert S == length(Q) == length(V) == length(Y)
    @assert D == length(Q[1]) == length(V[1]) == length(Y[1])

    # copy x to Y and compute Q = q + Δt Y
    for k in 1:D
        Y[i][k] = x[k]
        Q[i][k] = params.q[k] + params.Δt * Y[i][k]
    end

    # compute V = v(Q)
    tᵢ = params.t + params.Δt * params.tab.q.c[i]
    params.equ.v(tᵢ, Q[i], V[i])
end


"Compute stages of fully implicit Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDIRK{DT,TT,ET,D,S}, i::Int) where {ST,DT,TT,ET,D,S}

    cache = NonlinearFunctionCacheDIRK{ST}(D, S)

    quote
        compute_stages!(x, $cache.Q, $cache.V, $cache.Y, params, i)

        local y1::ST
        local y2::ST

        # compute b = - (Y-AV)
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * $cache.V[j][k]
                y2 += params.tab.q.â[i,j] * $cache.V[j][k]
            end
            b[k] = - $cache.Y[i][k] + (y1 + y2)
        end
    end
end


"Diagonally implicit Runge-Kutta integrator."
struct IntegratorDIRK{DT, TT, PT <: ParametersDIRK{DT,TT},
                              ST,# <: NonlinearSolver{DT},
                              IT <: InitialGuessODE{DT,TT}, N} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
end

function IntegratorDIRK(equation::ODE{DT,TT,FT,N}, tableau::TableauDIRK{TT}, Δt::TT) where {DT,TT,FT,N}
    D = equation.d
    M = equation.n
    S = tableau.q.s

    # create params
    params = ParametersDIRK(equation, tableau, Δt)

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
    iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorDIRK{DT, TT, typeof(params), typeof(solvers), typeof(iguess), N}(
                params, solvers, iguess)

end

equation(int::IntegratorDIRK) = int.params.equ
timestep(int::IntegratorDIRK) = int.params.Δt
has_initial_guess(int::IntegratorDIRK) = true

"""
Diagonally implicit Runge-Kutta integrator cache.
"""
mutable struct IntegratorCacheDIRK{DT,TT,D,S} <: ODEIntegratorCache{DT,D}
    n::Int
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}

    qₑᵣᵣ::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}

    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCacheDIRK{DT,TT,D,S}() where {DT,TT,D,S}
        q = zeros(TwicePrecision{DT}, D)
        q̅ = zeros(TwicePrecision{DT}, D)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        new(0, zero(TT), zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Y)
    end
end

function create_integrator_cache(int::IntegratorDIRK{DT,TT}) where {DT,TT}
    IntegratorCacheDIRK{DT, TT, ndims(int), int.params.tab.s}()
end


"Initialise initial guess"
function initialize!(int::IntegratorDIRK, cache::IntegratorCacheDIRK)
    # initialise initial guess
    cache.t̅ = cache.t - timestep(int)

    int.params.equ.v(cache.t, cache.q, cache.v)

    initialize!(int.iguess, cache.t, cache.q, cache.v, cache.t̅, cache.q̅, cache.v̅)
end


"Compute initial guess for internal stages."
function initial_guess!(int::IntegratorDIRK, cache::IntegratorCacheDIRK)
    for i in 1:int.params.tab.q.s
        evaluate!(int.iguess, cache.q, cache.v, cache.q̅, cache.v̅, cache.q̃, cache.ṽ, int.params.tab.q.c[i])
        for k in eachindex(cache.V[i], cache.ṽ)
            cache.V[i][k] = cache.ṽ[k]
        end
    end
    for i in 1:int.params.tab.q.s
        for k in 1:ndims(int)
            int.solver[i].x[k] = 0
            for j in 1:int.params.tab.q.s
                int.solver[i].x[k] += int.params.tab.q.a[i,j] * cache.V[j][k]
            end
        end
    end
end


"Integrate ODE with diagonally implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorDIRK{DT,TT}, cache::IntegratorCacheDIRK{DT,TT}) where {DT,TT,N}
    int.params.t  = cache.t
    int.params.q .= cache.q

    # compute initial guess
    initial_guess!(int, cache)

    # reset cache
    reset!(cache, timestep(int))

    # consecutively solve for all stages
    for i in 1:nstages(int)
        # call nonlinear solver
        solve!(int.solver[i])

        # print solver status
        print_solver_status(int.solver[i].status, int.solver[i].params)

        # check if solution contains NaNs or error bounds are violated
        check_solver_status(int.solver[i].status, int.solver[i].params)

        # compute vector field at internal stages
        compute_stages!(int.solver[i].x, cache.Q, cache.V, cache.Y, int.params, i)
    end

    # compute final update
    update_solution!(cache.q, cache.qₑᵣᵣ, cache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)

    # update vector field for initial guess
    update!(int.iguess, cache.t, cache.q, cache.v)

    # take care of periodic solutions
    cut_periodic_solution!(cache, int.params.equ.periodicity)
end

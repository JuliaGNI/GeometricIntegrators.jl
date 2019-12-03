
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
struct IntegratorDIRK{DT, TT, PT <: ParametersDIRK{DT,TT},
                              ST,# <: NonlinearSolver{DT},
                              IT <: InitialGuessODE{DT,TT}, N, D, S} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheDIRK{DT,D,S}
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

    # create cache
    cache = IntegratorCacheDIRK{DT,D,S}()

    # create integrator
    IntegratorDIRK{DT, TT, typeof(params), typeof(solvers), typeof(iguess), N, D, S}(
                params, solvers, iguess, cache)

end

has_initial_guess(int::IntegratorDIRK) = true


"Initialise initial guess"
function initialize!(int::IntegratorDIRK, cache::IntegratorCacheDIRK)
    # initialise initial guess
    cache.t̅ = cache.t - timestep(int)

    int.params.equ.v(cache.t, cache.q, cache.v)

    initialize!(int.iguess, cache.t, cache.q, cache.v, cache.t̅, cache.q̅, cache.v̅)
end


"Compute initial guess for internal stages."
function initial_guess!(int::IntegratorDIRK, sol::AtomisticSolutionODE)
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

    cache = IntegratorCacheDIRK{ST,D,S}()

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


"Integrate ODE with diagonally implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorDIRK{DT,TT}, sol::AtomisticSolutionODE{DT,TT}) where {DT,TT,N}
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
        compute_stages!(int.solver[i].x, int.cache.Q, int.cache.V, int.cache.Y, int.params, i)
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)

    # update vector field for initial guess
    update!(int.iguess, sol.t, sol.q, sol.v)
end

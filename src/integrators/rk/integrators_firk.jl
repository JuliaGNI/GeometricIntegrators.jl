
"Holds the tableau of a fully implicit Runge-Kutta method."
struct TableauFIRK{T} <: AbstractTableauIRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauFIRK{T}(q) where {T}
        if (q.s > 1 && istrilstrict(q.a)) || (q.s==1 && q.a[1,1] == 0)
            @warn "Initializing TableauFIRK with explicit tableau $(q.name).\nYou might want to use TableauERK instead."
        elseif q.s > 1 && istril(q.a)
            @warn "Initializing TableauFIRK with diagonally implicit tableau $(q.name).\nYou might want to use TableauDIRK instead."
        end

        new(q.name, q.o, q.s, q)
    end
end

function TableauFIRK(q::CoefficientsRK{T}) where {T}
    TableauFIRK{T}(q)
end

function TableauFIRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauFIRK{T}(CoefficientsRK(name, order, a, b, c))
end

# TODO function readTableauFIRKFromFile(dir::AbstractString, name::AbstractString)


"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct ParametersFIRK{DT, TT, ET <: ODE{DT,TT}, D, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauFIRK{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
end

function ParametersFIRK(equ::ET, tab::TableauFIRK{TT}, Δt::TT) where {DT, TT, ET <: ODE{DT,TT}}
    ParametersFIRK{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, 0, zeros(DT, equ.d))
end


"""
Fully implicit Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of solution
* `ṽ`: initial guess of vector field
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Y`: vector field of internal stages
"""
struct IntegratorCacheFIRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function IntegratorCacheFIRK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D), Q, V, Y)
    end
end

"Fully implicit Runge-Kutta integrator."
struct IntegratorFIRK{DT, TT, PT <: ParametersFIRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessODE{DT,TT}, N, D, S} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheFIRK{DT,D,S}
end

function IntegratorFIRK(equation::ODE{DT,TT,FT,N}, tableau::TableauFIRK{TT}, Δt::TT) where {DT,TT,FT,N}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersFIRK(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*S, params)

    # create initial guess
    iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheFIRK{DT,D,S}()

    # create integrator
    IntegratorFIRK{DT, TT, typeof(params), typeof(solver), typeof(iguess), N, D, S}(
                params, solver, iguess, cache)
end


function initialize!(int::IntegratorFIRK, sol::AtomisticSolutionODE)
    sol.t̅ = sol.t - timestep(int)

    int.params.equ.v(sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̅, sol.q̅, sol.v̅)
end


function update_params!(int::IntegratorFIRK, sol::AtomisticSolutionODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
end


function initial_guess!(int::IntegratorFIRK, sol::AtomisticSolutionODE)
    local offset::Int

    # compute initial guess for internal stages
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.v, sol.q̅, sol.v̅, int.cache.q̃, int.cache.ṽ, tableau(int).q.c[i])
        for k in eachindex(int.cache.V[i], int.cache.ṽ)
            int.cache.V[i][k] = int.cache.ṽ[k]
        end
    end
    for i in eachstage(int)
        offset = ndims(int)*(i-1)
        for k in eachdim(int)
            int.solver.x[offset+k] = 0
            for j in eachstage(int)
                int.solver.x[offset+k] += tableau(int).q.a[i,j] * int.cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                                    params::ParametersFIRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    local tᵢ::TT

    @assert S == length(Q) == length(V) == length(Y)

    # copy x to Y and compute Q = q + Δt Y
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i]) == length(Y[i])
        for k in 1:D
            Y[i][k] = x[D*(i-1)+k]
            Q[i][k] = params.q[k] + params.Δt * Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in 1:S
        tᵢ = params.t + params.Δt * params.tab.q.c[i]
        params.equ.v(tᵢ, Q[i], V[i])
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersFIRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    cache = IntegratorCacheFIRK{ST,D,S}()

    quote
        compute_stages!(x, $cache.Q, $cache.V, $cache.Y, params)

        local y1::ST
        local y2::ST

        # compute b = - (Y-AV)
        for i in 1:S
            for k in 1:D
                y1 = 0
                y2 = 0
                for j in 1:S
                    y1 += params.tab.q.a[i,j] * $cache.V[j][k]
                    y2 += params.tab.q.â[i,j] * $cache.V[j][k]
                end
                b[D*(i-1)+k] = - $cache.Y[i][k] + (y1 + y2)
            end
        end
    end
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFIRK{DT,TT}, sol::AtomisticSolutionODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.Y, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.v)
end

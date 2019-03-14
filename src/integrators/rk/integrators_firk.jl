
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

struct NonlinearFunctionCacheFIRK{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    v::Vector{DT}
    y::Vector{DT}

    function NonlinearFunctionCacheFIRK{DT}(D,S) where {DT}

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
                                    params::ParametersFIRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}

    local tᵢ::TT

    @assert S == length(Q) == length(V) == length(Y)
    @assert D == length(Q[1]) == length(V[1]) == length(Y[1])

    # copy x to Y and compute Q = q + Δt Y
    for i in 1:S
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

    cache = NonlinearFunctionCacheFIRK{ST}(D, S)

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


"Fully implicit Runge-Kutta integrator."
struct IntegratorFIRK{DT, TT, PT <: ParametersFIRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessODE{DT,TT}, N} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    fcache::NonlinearFunctionCacheFIRK{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
end

function IntegratorFIRK(equation::ODE{DT,TT,FT,N}, tableau::TableauFIRK{TT}, Δt::TT) where {DT,TT,FT,N}
    D = equation.d
    M = equation.n
    S = tableau.q.s

    # create params
    params = ParametersFIRK(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*S, params)

    # create initial guess
    iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    fcache = NonlinearFunctionCacheFIRK{DT}(D, S)

    # create solution vectors
    q = create_solution_vector(DT, D, M)

    # create integrator
    IntegratorFIRK{DT, TT, typeof(params), typeof(solver), typeof(iguess), N}(
                params, solver, iguess, fcache, q)
end


function initialize!(int::IntegratorFIRK{DT,TT}, sol::SolutionODE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m])
end

function initial_guess!(int::IntegratorFIRK, m::Int)
    local offset::Int

    # compute initial guess for internal stages
    for i in 1:int.params.tab.q.s
        evaluate!(int.iguess, m, int.fcache.y, int.fcache.v, int.params.tab.q.c[i])
        for k in eachindex(int.fcache.V[i], int.fcache.v)
            int.fcache.V[i][k] = int.fcache.v[k]
        end
    end
    for i in 1:int.params.tab.q.s
        offset = dims(int)*(i-1)
        for k in 1:dims(int)
            int.solver.x[offset+k] = 0
            for j in 1:int.params.tab.q.s
                int.solver.x[offset+k] += int.params.tab.q.a[i,j] * int.fcache.V[j][k]
            end
        end
    end
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFIRK{DT,TT}, sol::SolutionODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,N}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.params.Δt
    int.params.q .= int.q[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, int.fcache.Q, int.fcache.V, int.fcache.Y, int.params)

    # compute final update
    update_solution!(int.q[m], int.fcache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)

    # copy solution to initial guess
    update!(int.iguess, m, sol.t[0] + n*int.params.Δt, int.q[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], n, m)
end

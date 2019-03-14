@doc raw"""
`TableauIPRK`: Tableau of an Implicit Partitioned Runge-Kutta method
```math
\begin{align*}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} ,
\end{align*}
```
usually satisfying the symplecticity conditions
```math
\begin{align*}
b_{i} \bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\bar{b}_i &= b_i .
\end{align*}
```
"""
struct TableauIPRK{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}
    p::CoefficientsRK{T}

    function TableauIPRK{T}(name, o, q, p) where {T}
        @assert q.s==p.s
        new(name, o, q.s, q, p)
    end
end

function TableauIPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}) where {T}
    TableauIPRK{T}(name, order, q, p)
end

function TableauIPRK(name::Symbol, order::Int, q::CoefficientsRK{T}) where {T}
    TableauIPRK{T}(name, order, q, q)
end

# TODO function readTableauIPRKFromFile(dir::AbstractString, name::AbstractString)


"Parameters for right-hand side function of implicit partitioned Runge-Kutta methods."
mutable struct ParametersIPRK{DT, TT, ET <: PODE{DT,TT}, D, S} <: Parameters{DT,TT}
    equ::ET
    tab::TableauIPRK{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersIPRK(equ::ET, tab::TableauIPRK{TT}, Δt::TT) where {DT, TT, ET <: PODE{DT,TT}}
    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    ParametersIPRK{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, 0, q, p)
end


struct NonlinearFunctionCacheIPRK{ST}
    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}
    Z::Vector{Vector{ST}}

    v::Vector{ST}
    f::Vector{ST}
    y::Vector{ST}
    z::Vector{ST}

    function NonlinearFunctionCacheIPRK{ST}(D,S) where {ST}
        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        # create update vectors
        v = zeros(ST,D)
        f = zeros(ST,D)
        y = zeros(ST,D)
        z = zeros(ST,D)

        new(Q, V, P, F, Y, Z, v, f, y, z)
    end
end

"Compute stages of implicit partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersIPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    cache = NonlinearFunctionCacheIPRK{ST}(D, S)

    quote
        compute_stages!(x, $cache.Q, $cache.V, $cache.P, $cache.F, $cache.Y, $cache.Z, params)

        # compute b = - [(Y-AV), (Z-AF)]
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+1] = - $cache.Y[i][k]
                b[2*(D*(i-1)+k-1)+2] = - $cache.Z[i][k]
                for j in 1:S
                    b[2*(D*(i-1)+k-1)+1] += params.tab.q.a[i,j] * $cache.V[j][k]
                    b[2*(D*(i-1)+k-1)+2] += params.tab.p.a[i,j] * $cache.F[j][k]
                end
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}},
                                        P::Vector{Vector{ST}}, F::Vector{Vector{ST}},
                                        Y::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        params::ParametersIPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    local tqᵢ::TT
    local tpᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy y to Y and Z
            Y[i][k] = x[2*(D*(i-1)+k-1)+1]
            Z[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            Q[i][k] = params.q[k] + params.Δt * Y[i][k]
            P[i][k] = params.p[k] + params.Δt * Z[i][k]
        end

        # compute time of internal stage
        tqᵢ = params.t + params.Δt * params.tab.q.c[i]
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]

        # compute v(Q,P) and f(Q,P)
        params.equ.v(tqᵢ, Q[i], P[i], V[i])
        params.equ.f(tpᵢ, Q[i], P[i], F[i])
    end
end


"Implicit partitioned Runge-Kutta integrator."
struct IntegratorIPRK{DT, TT, PT <: ParametersIPRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessPODE{DT,TT}} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    fcache::NonlinearFunctionCacheIPRK{DT}

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}
end

function IntegratorIPRK(equation::PODE{DT,TT,VT,FT}, tableau::TableauIPRK{TT}, Δt::TT) where {DT,TT,VT,FT}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersIPRK(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, 2*D*S, params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache for internal stage vectors and update vectors
    fcache = NonlinearFunctionCacheIPRK{DT}(D,S)

    # create solution vectors
    q = create_solution_vector(DT, D, M)
    p = create_solution_vector(DT, D, M)

    # create integrator
    IntegratorIPRK{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(
                params, solver, iguess, fcache, q, p)
end

equation(integrator::IntegratorIPRK) = integrator.params.equ
timestep(integrator::IntegratorIPRK) = integrator.params.Δt


function initialize!(int::IntegratorIPRK, sol::SolutionPODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m], int.p[m])
end

function initial_guess!(int::IntegratorIPRK, m::Int)
    for i in 1:int.params.tab.q.s
        evaluate!(int.iguess, m, int.fcache.y, int.fcache.z, int.fcache.v, int.fcache.f, int.params.tab.q.c[i], int.params.tab.p.c[i])
        for k in 1:int.params.equ.d
            int.fcache.V[i][k] = int.fcache.v[k]
            int.fcache.F[i][k] = int.fcache.f[k]
        end
    end
    for i in 1:int.params.tab.q.s
        for k in 1:int.params.equ.d
            int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+1] = 0
            int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+2] = 0
            for j in 1:int.params.tab.q.s
                int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+1] += int.params.tab.q.a[i,j] * int.fcache.V[j][k]
                int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+2] += int.params.tab.p.a[i,j] * int.fcache.F[j][k]
            end
        end
    end
end

"Integrate ODE with implicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorIPRK{DT,TT}, sol::SolutionPODE{DT,TT}, m::Int, n::Int) where {DT,TT}
    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.params.Δt

    # copy previous solution
    int.params.q .= int.q[m]
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.fcache.Q, int.fcache.V, int.fcache.P, int.fcache.F, int.fcache.Y, int.fcache.Z, int.params)

    # compute final update
    update_solution!(int.q[m], int.fcache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(int.p[m], int.fcache.F, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.params.equ.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end

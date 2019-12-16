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


"""
Implicit partitioned Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of q
* `p̃`: initial guess of p
* `ṽ`: initial guess of v
* `f̃`: initial guess of f
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of q
* `P`: internal stages of p
* `V`: internal stages of v
* `F`: internal stages of f
* `Y`: vector field of internal stages of q
* `Z`: vector field of internal stages of p
"""
struct IntegratorCacheIPRK{ST,D,S} <: PODEIntegratorCache{ST,D}
    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    s̃::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}
    Z::Vector{Vector{ST}}

    function IntegratorCacheIPRK{ST,D,S}() where {ST,D,S}
        # create temporary vectors
        q̃ = zeros(ST,D)
        p̃ = zeros(ST,D)
        ṽ = zeros(ST,D)
        f̃ = zeros(ST,D)
        s̃ = zeros(ST,D)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        new(q̃, p̃, ṽ, f̃, s̃, Q, P, V, F, Y, Z)
    end
end


"Implicit partitioned Runge-Kutta integrator."
struct IntegratorIPRK{DT, TT, PT <: ParametersIPRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessPODE{DT,TT}, D, S} <: IntegratorPRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheIPRK{DT,D,S}
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

    # create cache
    cache = IntegratorCacheIPRK{DT,D,S}()

    # create integrator
    IntegratorIPRK{DT, TT, typeof(params), typeof(solver), typeof(iguess), D, S}(params, solver, iguess, cache)
end


function update_params!(int::IntegratorIPRK, sol::AtomisticSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.p .= sol.p
end


"Compute stages of implicit partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersIPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    cache = IntegratorCacheIPRK{ST,D,S}()

    quote
        compute_stages!(x, $cache, params)

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


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheIPRK{ST,D,S},
                                        params::ParametersIPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    local tqᵢ::TT
    local tpᵢ::TT

    for i in 1:S
        for k in 1:D
            # copy y to Y and Z
            cache.Y[i][k] = x[2*(D*(i-1)+k-1)+1]
            cache.Z[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            cache.Q[i][k] = params.q[k] + params.Δt * cache.Y[i][k]
            cache.P[i][k] = params.p[k] + params.Δt * cache.Z[i][k]
        end

        # compute time of internal stage
        tqᵢ = params.t + params.Δt * params.tab.q.c[i]
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]

        # compute v(Q,P) and f(Q,P)
        params.equ.v(tqᵢ, cache.Q[i], cache.P[i], cache.V[i])
        params.equ.f(tpᵢ, cache.Q[i], cache.P[i], cache.F[i])
    end
end


function initialize!(int::IntegratorIPRK, sol::AtomisticSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    int.params.equ.v(sol.t, sol.q, sol.p, sol.v)
    int.params.equ.f(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end


function initial_guess!(int::IntegratorIPRK, sol::AtomisticSolutionPODE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q.c[i],
                              tableau(int).p.c[i])

        for k in eachdim(int)
            int.cache.V[i][k] = int.cache.ṽ[k]
            int.cache.F[i][k] = int.cache.f̃[k]
        end
    end
    for i in eachstage(int)
        for k in eachdim(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] = 0
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] = 0
            for j in eachstage(int)
                int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] += tableau(int).q.a[i,j] * int.cache.V[j][k]
                int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] += tableau(int).p.a[i,j] * int.cache.F[j][k]
            end
        end
    end
end


"Integrate ODE with implicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorIPRK{DT,TT}, sol::AtomisticSolutionPODE{DT,TT}) where {DT,TT}
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

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

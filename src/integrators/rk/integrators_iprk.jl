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


"Implicit partitioned Runge-Kutta integrator."
struct IntegratorIPRK{DT, TT, PT <: ParametersIPRK{DT,TT},
                              ST <: NonlinearSolver{DT},
                              IT <: InitialGuessPODE{DT,TT}} <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
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

    # create integrator
    IntegratorIPRK{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(params, solver, iguess)
end

equation(integrator::IntegratorIPRK) = integrator.params.equ
timestep(integrator::IntegratorIPRK) = integrator.params.Δt
has_initial_guess(int::IntegratorIPRK) = true


"""
Implicit partitioned Runge-Kutta integrator cache.

### Fields

* `n`: time step number
* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `p`: current solution of p
* `p̅`: previous solution of p
* `v`: vector field of q
* `v̅`: vector field of q̅
* `f`: vector field of p
* `f̅`: vector field of p̅
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
mutable struct IntegratorCacheIPRK{ST,TT,D,S} <: ODEIntegratorCache{ST,D,S}
    n::Int
    t::TT
    t̅::TT

    q::Vector{TwicePrecision{ST}}
    q̅::Vector{TwicePrecision{ST}}
    p::Vector{TwicePrecision{ST}}
    p̅::Vector{TwicePrecision{ST}}

    v::Vector{ST}
    v̅::Vector{ST}
    f::Vector{ST}
    f̅::Vector{ST}

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

    function IntegratorCacheIPRK{ST,TT,D,S}() where {ST,TT,D,S}
        q = zeros(TwicePrecision{ST}, D)
        q̅ = zeros(TwicePrecision{ST}, D)
        p = zeros(TwicePrecision{ST}, D)
        p̅ = zeros(TwicePrecision{ST}, D)

        # create update vectors
        v = zeros(ST,D)
        v̅ = zeros(ST,D)
        f = zeros(ST,D)
        f̅ = zeros(ST,D)

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

        new(0, zero(TT), zero(TT), q, q̅, p, p̅, v, v̅, f, f̅, q̃, ṽ, p̃, f̃, s̃, Q, P, V, F, Y, Z)
    end
end

function create_integrator_cache(int::IntegratorIPRK{DT,TT}) where {DT,TT}
    IntegratorCacheIPRK{DT, TT, ndims(equation(int)), int.params.tab.s}()
end

function CommonFunctions.reset!(cache::IntegratorCacheIPRK{DT,TT}, Δt::TT) where {DT,TT}
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.p̅ .= cache.p
    cache.t += Δt
    cache.n += 1
end

function cut_periodic_solution!(cache::IntegratorCacheIPRK, periodicity::Vector)
    cut_periodic_solution!(cache.q, periodicity, cache.s̃)
    cache.q .+= cache.s̃
    cache.q̅ .+= cache.s̃
end

function CommonFunctions.get_solution(cache::IntegratorCacheIPRK)
    (cache.t, cache.q, cache.p)
end

function CommonFunctions.set_solution!(cache::IntegratorCacheIPRK, sol, n=0)
    t, q, p = sol
    cache.n  = n
    cache.t  = t
    cache.q .= q
    cache.p .= p
    cache.v .= 0
    cache.f .= 0
end


"Compute stages of implicit partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersIPRK{DT,TT,ET,D,S}) where {ST,DT,TT,ET,D,S}
    cache = IntegratorCacheIPRK{ST,TT,D,S}()

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


function compute_stages!(x::Vector{ST}, cache::IntegratorCacheIPRK{ST,TT,D,S},
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


function initial_guess!(int::IntegratorIPRK, cache::IntegratorCacheIPRK)
    for i in 1:int.params.tab.s
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              int.params.tab.q.c[i],
                              int.params.tab.p.c[i])

        for k in 1:int.params.equ.d
            cache.V[i][k] = cache.ṽ[k]
            cache.F[i][k] = cache.f̃[k]
        end
    end
    for i in 1:int.params.tab.s
        for k in 1:int.params.equ.d
            int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+1] = 0
            int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+2] = 0
            for j in 1:int.params.tab.s
                int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+1] += int.params.tab.q.a[i,j] * cache.V[j][k]
                int.solver.x[2*(int.params.equ.d*(i-1)+k-1)+2] += int.params.tab.p.a[i,j] * cache.F[j][k]
            end
        end
    end
end

"Integrate ODE with implicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorIPRK{DT,TT}, cache::IntegratorCacheIPRK{DT,TT}) where {DT,TT}
    # set time for nonlinear solver and copy previous solution
    int.params.t  = cache.t
    int.params.q .= cache.q
    int.params.p .= cache.p

    # compute initial guess
    initial_guess!(int, cache)

    # reset cache
    reset!(cache, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, cache.n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, cache.n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # compute final update
    update_solution!(cache.q, cache.V, int.params.tab.q.b, int.params.tab.q.b̂, int.params.Δt)
    update_solution!(cache.p, cache.F, int.params.tab.p.b, int.params.tab.p.b̂, int.params.Δt)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, int.params.equ.periodicity)
end

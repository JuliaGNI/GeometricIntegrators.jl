
"Parameters for right-hand side function of implicit partitioned Runge-Kutta methods."
mutable struct ParametersPRKimplicit{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::PartitionedTableau{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    p::Vector{DT}

    function ParametersPRKimplicit{DT,D}(equs::ET, tab::PartitionedTableau{TT}, Δt::TT) where {D, DT, TT, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end


@doc raw"""
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
struct IntegratorCachePRKimplicit{ST,D,S} <: PODEIntegratorCache{ST,D}
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

    function IntegratorCachePRKimplicit{ST,D,S}() where {ST,D,S}
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

function IntegratorCache(params::ParametersPRKimplicit{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCachePRKimplicit{DT,D,S}(; kwargs...)
end

function IntegratorCache{ST}(params::ParametersPRKimplicit{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCachePRKimplicit{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersPRKimplicit{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCachePRKimplicit{ST,D,S}


@doc raw"""
Implicit partitioned Runge-Kutta integrator solving the system
```math
\begin{aligned}
P_{n,i} &= \vartheta (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} ,
\end{aligned}
```
Usually we are interested in Lagrangian systems, where
```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) ,
\end{aligned}
```
and tableaus satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```
"""
struct IntegratorPRKimplicit{DT, TT, D, S, PT <: ParametersPRKimplicit{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessIODE{TT}} <: AbstractIntegratorIRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorPRKimplicit(params::ParametersPRKimplicit{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorPRKimplicit{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersPRKimplicit{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorPRKimplicit(params, solver, iguess, caches)
    end

    function IntegratorPRKimplicit{DT,D}(v::Function, f::Function, tableau::PartitionedTableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorPRKimplicit{DT,D}(NamedTuple{(:v,:f)}((v,f)), tableau, Δt; kwargs...)
    end

    function IntegratorPRKimplicit{DT,D}(v::Function, f::Function, h::Function, tableau::PartitionedTableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorPRKimplicit{DT,D}(NamedTuple{(:v,:f,:h)}((v,f,h)), tableau, Δt; kwargs...)
    end

    function IntegratorPRKimplicit(problem::Union{IODEProblem{DT}, LODEProblem{DT}}, tableau::PartitionedTableau{TT}; kwargs...) where {DT,TT}
        IntegratorPRKimplicit{DT, ndims(problem)}(functions(problem), tableau, timestep(problem); kwargs...)
    end
end


@inline Base.ndims(::IntegratorPRKimplicit{DT,TT,D,S}) where {DT,TT,D,S} = D


function initialize!(int::IntegratorPRKimplicit, sol::AtomicSolutionPODE)
    sol.t̄ = sol.t - timestep(int)

    equations(int)[:v̄](sol.t, sol.q, sol.v)
    equations(int)[:f̄](sol.t, sol.q, sol.v, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end


function update_params!(int::IntegratorPRKimplicit, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.p .= sol.p
end


function initial_guess!(int::IntegratorPRKimplicit{DT}, sol::AtomicSolutionPODE{DT},
                        cache::IntegratorCachePRKimplicit{DT}=int.caches[DT]) where {DT}

    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.Q[i], cache.P[i], cache.V[i], cache.F[i],
                              tableau(int).q.c[i], tableau(int).p.c[i])
    end
    for i in eachstage(int)
        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = 0
            for j in eachstage(int)
                int.solver.x[ndims(int)*(i-1)+k] += tableau(int).q.a[i,j] * cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                                        P::Vector{Vector{ST}}, F::Vector{Vector{ST}}, Z::Vector{Vector{ST}},
                                        params::ParametersPRKimplicit{DT,TT,D,S}) where {ST,DT,TT,D,S}
    local tqᵢ::TT
    local tpᵢ::TT
    local y1::ST
    local y2::ST

    for i in 1:S
        for k in 1:D
            # copy x to V 
            V[i][k] = x[D*(i-1)+k]

            # compute Q
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
            end
            Q[i][k] = params.q[k] + params.Δt * (y1 + y2)
        end

        # compute time of internal stage
        tqᵢ = params.t + params.Δt * params.tab.q.c[i]
        tpᵢ = params.t + params.Δt * params.tab.p.c[i]

        # compute ϑ(Q,V) and f(Q,V)
        params.equs[:ϑ](tqᵢ, Q[i], V[i], P[i])
        params.equs[:f](tpᵢ, Q[i], V[i], F[i])
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersPRKimplicit{DT,TT,D,S},
                          caches::CacheDict) where {ST,DT,TT,D,S}

    # temporary variables
    local z1::ST
    local z2::ST

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q, cache.V, cache.Y, cache.P, cache.F, cache.Z, params)

    # compute b = [P-p-AF]
    for i in 1:S
        for k in 1:D
            z1 = 0
            z2 = 0
            for j in 1:S
                z1 += params.tab.p.a[i,j] * cache.F[j][k]
                z2 += params.tab.p.â[i,j] * cache.F[j][k]
            end
            b[D*(i-1)+k] = cache.P[i][k] - params.p[k] - params.Δt * (z1 + z2)
        end
    end
end


function integrate_step!(int::IntegratorPRKimplicit{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCachePRKimplicit{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from atomic solution
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.Y,
                                  cache.P, cache.F, cache.Z, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

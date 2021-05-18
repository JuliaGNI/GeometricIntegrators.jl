
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct ParametersSRKimplicit{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::Tableau{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    θ::Vector{DT}

    function ParametersSRKimplicit{DT,D}(equs::ET, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end


@doc raw"""
Fully implicit Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of solution
* `ṽ`: initial guess of vector field
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Θ`: implicit function of internal stages
* `F`: vector field of implicit function
"""
struct IntegratorCacheSRKimplicit{DT,D,S} <: ODEIntegratorCache{DT,D}
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Θ::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function IntegratorCacheSRKimplicit{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D),
            Q, V, Θ, F)
    end
end

function IntegratorCache{ST}(params::ParametersSRKimplicit{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheSRKimplicit{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersSRKimplicit{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheSRKimplicit{ST,D,S}


@doc raw"""
Special implicit Runge-Kutta integrator for certain noncanonical
symplectic equations solving the system
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
struct IntegratorSRKimplicit{DT, TT, D, S, PT <: ParametersSRKimplicit{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorSRKimplicit(params::ParametersSRKimplicit{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorSRKimplicit{DT,D}(equations::NamedTuple, tableau::Tableau{TT}, Δt::TT; exact_jacobian=true) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersSRKimplicit{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessODE(get_config(:ig_extrapolation), equations[:v̄], Δt)

        # create integrator
        IntegratorSRKimplicit(params, solver, iguess, caches)
    end

    # function IntegratorSRKimplicit{DT,D}(v::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
    #     IntegratorSRKimplicit{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    # end

    # function IntegratorSRKimplicit{DT,D}(v::Function, h::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
    #     IntegratorSRKimplicit{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    # end

    function IntegratorSRKimplicit(equation::Union{IODE{DT}, LODE{DT}}, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorSRKimplicit{DT, ndims(equation)}(get_functions(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorSRKimplicit{DT,TT,D,S}) where {DT,TT,D,S} = D


Solutions.AtomicSolution(integrator::IntegratorSRKimplicit{DT,TT}) where {DT,TT} =
    AtomicSolutionPODE(DT, TT, ndims(integrator), get_internal_variables(integrator))


function initialize!(int::IntegratorSRKimplicit, sol::AtomicSolutionODE)
    sol.t̄ = sol.t - timestep(int)

    equations(int)[:v̄](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)
end


function update_params!(int::IntegratorSRKimplicit, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    equations(int)[:ϑ](sol.t, sol.q, sol.v, int.params.θ)
end


function initial_guess!(int::IntegratorSRKimplicit{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheSRKimplicit{DT}=int.caches[DT]) where {DT,TT}

    # compute initial guess for internal stages
    # for i in eachstage(int)
    #     evaluate!(int.iguess, sol.q̄, sol.v̄, sol.q, sol.v, cache.Q[i], cache.V[i], tableau(int).c[i])
    # end
    for i in eachstage(int)
        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = 0#cache.V[i][k]
        end
    end

    # compute initial guess for solution
    # evaluate!(int.iguess, sol.q̄, sol.v̄, sol.q, sol.v, cache.q, cache.v, one(TT))
    for k in eachdim(int)
        int.solver.x[ndims(int)*nstages(int)+k] = sol.q[k]#cache.q[k]
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, 
                                        Θ::Vector{Vector{ST}}, F::Vector{Vector{ST}}, 
                                        q::Vector{ST}, v::Vector{ST}, θ::Vector{ST},
                                        params::ParametersSRKimplicit{DT,TT,D,S}) where {ST,DT,TT,D,S}

    # temporary variables
    local y1::ST
    local y2::ST
    local tᵢ::TT

    # copy x to V
    for i in eachindex(V)
        for k in eachindex(V[i])
            V[i][k] = x[D*(i-1)+k]
        end
    end
    for k in eachindex(q)
        q[k] = x[D*S+k]
    end

    # compute Q = q + Δt A V, Θ = ϑ(Q), F = f(Q,V)
    for i in eachindex(Q,F,Θ)
        tᵢ = params.t + params.Δt * params.tab.c[i]
        for k in eachindex(Q[i],params.q)
            y1 = 0
            y2 = 0
            for j in eachindex(V)
                y1 += params.tab.a[i,j] * V[j][k]
                y2 += params.tab.â[i,j] * V[j][k]
            end
            Q[i][k] = params.q[k] + params.Δt * (y1 + y2)
        end
        params.equs[:ϑ](tᵢ, Q[i], V[i], Θ[i])
        params.equs[:f](tᵢ, Q[i], V[i], F[i])
    end

    # compute q̄ = q + Δt B V, Θ = ϑ(q̄)
    tᵢ = params.t + params.Δt
    params.equs[:ϑ](tᵢ, q, v, θ)
end


# Compute stages of fully implicit Runge-Kutta methods.
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersSRKimplicit{DT,TT,D,S},
                          caches::CacheDict) where {ST,DT,TT,D,S}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q, cache.V, cache.Θ, cache.F, cache.q, cache.v, cache.θ, params)

    # compute b
    for i in eachindex(cache.Θ)
        for k in eachindex(cache.Θ[i], params.θ)
            y1 = 0
            y2 = 0
            for j in eachindex(cache.F)
                y1 += params.tab.a[i,j] * cache.F[j][k]
                y2 += params.tab.â[i,j] * cache.F[j][k]
            end
            b[D*(i-1)+k] = cache.Θ[i][k] - params.θ[k] - params.Δt * (y1 + y2)
        end
    end
    for k in 1:div(D,2)
        y1 = 0
        y2 = 0
        for j in eachindex(cache.F)
            y1 += params.tab.b[j] * cache.F[j][k]
            y2 += params.tab.b̂[j] * cache.F[j][k]
        end
        b[D*S+k] = cache.θ[k] - params.θ[k] - params.Δt * (y1 + y2)
    end
    for k in 1:div(D,2)
        y1 = 0
        y2 = 0
        for j in eachindex(cache.V)
            y1 += params.tab.b[j] * cache.V[j][k]
            y2 += params.tab.b̂[j] * cache.V[j][k]
        end
        b[D*S+div(D,2)+k] = cache.q[k] - params.q[k] - params.Δt * (y1 + y2)
    end
end


function integrate_step!(int::IntegratorSRKimplicit{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheSRKimplicit{DT}=int.caches[DT]) where {DT,TT}

    # update nonlinear solver parameters from atomic solution
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.Θ, cache.F, cache.q, cache.v, cache.θ, int.params)

    # compute final update
    sol.q .= cache.q
    sol.p .= cache.θ

    # compute vector field for initial guess
    equations(int)[:v̄](sol.t, sol.q, sol.v)
    # update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end


"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRK = AbstractParametersVPRK{:vprk}

function IntegratorCache(params::ParametersVPRK{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(false; kwargs...)
end

function IntegratorCache{ST}(params::ParametersVPRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(false; kwargs...)
end


@doc raw"""
Variational Partitioned Runge-Kutta Integrator.

```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} - d_i \lambda , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} , \\
&&
0 &= \sum \limits_{i=1}^{s} d_i V_i , &&
\end{aligned}
```
satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```
"""
struct IntegratorVPRK{DT, TT, D, S, PT <: ParametersVPRK{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRK(params::ParametersVPRK{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRK{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, nullvec, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRK{DT,D}(equations, tableau, nullvec, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRK(params, solver, iguess, caches)
    end

    function IntegratorVPRK{DT,D}(ϑ::Function, f::Function, g::Function, v::Function,
                    tableau::PartitionedTableau{TT}, nullvec, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorVPRK{DT,D}(NamedTuple{(:ϑ,:f,:g,:v)}((ϑ, f, g, v)), tableau, nullvec, Δt; kwargs...)
    end

    function IntegratorVPRK{DT,D}(ϑ::Function, f::Function, g::Function, v::Function, h::Function,
                    tableau::PartitionedTableau{TT}, nullvec, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorVPRK{DT,D}(NamedTuple{(:ϑ,:f,:g,:v,:h)}((ϑ, f, g, v, h)), tableau, nullvec, Δt; kwargs...)
    end

    function IntegratorVPRK(problem::Union{IODEProblem{DT},LODEProblem{DT}}, tableau, nullvec; kwargs...) where {DT}
        IntegratorVPRK{DT, ndims(problem)}(functions(problem), tableau, nullvec, timestep(problem); kwargs...)
    end
end


IntegratorVPRKpNone = IntegratorVPRK


function Base.show(io::IO, int::IntegratorVPRK)
    print(io, "\nVariational Partitioned Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(int.params.Δt)\n")
    print(io, "   Tableau:  $(description(int.params.tab))\n")
    print(io, "   $(string(int.params.tab.q))")
    print(io, "   $(string(int.params.tab.p))")
    # print(io, reference(int.params.tab))
end


function Integrators.get_internal_variables(int::IntegratorVPRK{DT,TT,D,S}) where {DT, TT, D, S}
    Q = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)

    solver = get_solver_status(int.solver)

    (Q=Q, P=P, V=V, F=F, solver=solver)
end


function initial_guess!(int::IntegratorVPRK{DT}, sol::SolutionStepPODE{DT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRK{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.Q, cache.V, cache.P, cache.F, params)
    compute_rhs_vprk!(b, cache.P, cache.F, params)
    compute_rhs_vprk_correction!(b, cache.V, params)
end


function Integrators.integrate_step!(int::IntegratorVPRK{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset solution
    reset!(sol)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.P, cache.F, int.params)

    # compute final update
    update_solution!(int, sol, cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)

    # copy internal stage variables
    sol.internal[:Q] .= cache.Q
    sol.internal[:P] .= cache.P
    sol.internal[:V] .= cache.V
    sol.internal[:F] .= cache.F

    # copy solver status
    get_solver_status!(int.solver, sol.internal[:solver])
end

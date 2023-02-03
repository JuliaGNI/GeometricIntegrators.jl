
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpInternal = AbstractParametersVPRK{:vprk_pinternal}

function IntegratorCache(params::ParametersVPRKpInternal{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(D*(S+1), true; kwargs...)
end

function IntegratorCache{ST}(params::ParametersVPRKpInternal{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(D*(S+1), true; kwargs...)
end


"Variational Partitioned Runge-Kutta Integrator with Projection on Internal Stages."
struct IntegratorVPRKpInternal{DT, TT, D, S,
                PT <: ParametersVPRKpInternal{DT,TT},
                ST <: NonlinearSolver,
                IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    solver::ST
    iguess::IT
    caches::OldCacheDict{PT}

    function IntegratorVPRKpInternal(params::ParametersVPRKpInternal{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpInternal{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, nullvec, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        R = TT[1, tableau.R∞]
        params = ParametersVPRKpInternal{DT,D}(equations, tableau, nullvec, Δt, NamedTuple{(:R,)}((R,)))

        # create cache dict
        caches = OldCacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpInternal(params, solver, iguess, caches)
    end

    function IntegratorVPRKpInternal(problem::Union{IODEProblem{DT},LODEProblem{DT}}, tableau, nullvec; kwargs...) where {DT}
        IntegratorVPRKpInternal{DT, ndims(problem)}(functions(problem), tableau, nullvec, timestep(problem); kwargs...)
    end
end


function Base.show(io::IO, int::IntegratorVPRKpInternal)
    print(io, "\nVariational Partitioned Runge-Kutta Integrator with Projection on Internal Stages and:\n")
    print(io, "   Timestep: $(int.params.Δt)\n")
    print(io, "   Tableau:  $(description(int.params.tab))\n")
    print(io, "   $(string(int.params.tab.q))")
    print(io, "   $(string(int.params.tab.p))")
    # print(io, reference(int.params.tab))
end


function Integrators.get_internal_variables(int::IntegratorVPRKpInternal{DT,TT,D,S}) where {DT, TT, D, S}
    Q = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)
    λ = zeros(DT,D)

    solver = get_solver_status(int.solver)

    (Q=Q, P=P, V=V, F=F, λ=λ, solver=solver)
end


function initial_guess!(int::IntegratorVPRKpInternal{DT,TT}, sol::SolutionStepPODE{DT,TT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}

    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄[2], sol.p̄[2], sol.v̄[2], sol.f̄[2],
                              sol.q̄[1], sol.p̄[1], sol.v̄[1], sol.f̄[1],
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            cache.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end

    for k in eachdim(int)
        cache.x[ndims(int)*nstages(int)+k] = 0
    end
end


function compute_projection_vprk!(x::Vector{ST},
                q::SolutionVector{ST}, p::SolutionVector{ST}, v::SolutionVector{ST}, λ::SolutionVector{ST},
                Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpInternal{DT,TT,D,S}) where {ST,DT,TT,D,S}

    # create temporary variables
    local t₀::TT = params.t̄
    local t₁::TT = params.t̄ + params.Δt
    local tₘ::TT
    local y1::ST
    local y2::ST
    local g = zeros(ST,D)


    # copy x to λ and q̄
    for k in 1:D
        λ[k] = x[D*S+k]
    end

    # compute U=λ
    U[1] .= λ
    U[2] .= λ

    # compute Q
    for i in 1:S
        @assert D == length(Q[i]) == length(V[i])
        for k in 1:D
            y1 = 0
            y2 = 0
            for j in 1:S
                y1 += params.tab.q.a[i,j] * V[j][k]
                y2 += params.tab.q.â[i,j] * V[j][k]
            end
            Q[i][k] = params.q̄[k] + params.Δt * (y1 + y2) + params.Δt * params.pparams[:R][1] * U[1][k]
        end
    end

    # compute q
    for k in 1:D
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += params.tab.q.b[j] * V[j][k]
            y2 += params.tab.q.b̂[j] * V[j][k]
        end
        q[k] = params.q̄[k] + params.Δt * (y1 + y2) + params.Δt * (params.pparams[:R][1] * U[1][k] + params.pparams[:R][2] * U[2][k])
    end

    G[1] .= 0
    G[2] .= 0
    for j in 1:S
        tₘ = params.t̄ + params.Δt * params.tab.q.c[j]
        params.equ.g(g, tₘ, Q[j], V[j], λ)
        G[1] .+= params.tab.q.b[j] * g
        G[2] .+= params.tab.q.b[j] * g
    end

    # compute p=ϑ(q)
    params.equ.ϑ(p, t₁, q, v)
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpInternal{DT,TT,D,S},
                caches::OldCacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.q̃, cache.p̃, cache.ṽ, cache.λ, cache.Q, cache.V, cache.U, cache.P, cache.F, cache.G, params)

    # compute b = - [P-AF-U]
    compute_rhs_vprk!(b, cache.P, cache.F, cache.G, params)

    # compute b = - [p-bF-G]
    compute_rhs_vprk_projection_p!(b, cache.p̃, cache.F, cache.G, D*(S+0), params)

    compute_rhs_vprk_correction!(b, cache.V, params)
end


function Integrators.integrate_step!(int::IntegratorVPRKpInternal{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # reset cache
    reset!(sol, timestep(int))

    # compute initial guess
    initial_guess!(int, sol, cache)

    # call nonlinear solver
    solve!(cache.x, int.solver)

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages and projection vector fields
    compute_stages!(cache.x,
                    cache.q̃, cache.p̃, cache.ṽ, cache.λ,
                    cache.Q, cache.V, cache.U,
                    cache.P, cache.F, cache.G, int.params)

    # compute unprojected solution
    update_solution!(int, sol, cache)

    # add projection to solution
    project_solution!(int, sol, int.params.pparams[:R], cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)

    # copy internal stage variables
    # sol.internal.Q .= cache.Q
    # sol.internal.P .= cache.P
    # sol.internal.V .= cache.V
    # sol.internal.F .= cache.F
    # sol.internal.λ .= cache.λ

    # copy solver status
    # get_solver_status!(int.solver, sol.internal[:solver])
end

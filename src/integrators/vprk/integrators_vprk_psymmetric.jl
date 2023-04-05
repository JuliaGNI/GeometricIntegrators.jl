
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpSymmetric = AbstractParametersVPRK{:vprk_psymmetric}

function IntegratorCache(params::ParametersVPRKpSymmetric{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(D*(S+1), true; kwargs...)
end

function IntegratorCache{ST}(params::ParametersVPRKpSymmetric{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(D*(S+1), true; kwargs...)
end


"Variational Partitioned Runge-Kutta Integrator with Symmetric Projection."
struct IntegratorVPRKpSymmetric{DT, TT, D, S,
                PT <: ParametersVPRKpSymmetric{DT,TT},
                ST <: NonlinearSolver,
                IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    solver::ST
    iguess::IT
    caches::OldCacheDict{PT}

    function IntegratorVPRKpSymmetric(params::ParametersVPRKpSymmetric{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpSymmetric{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, nullvec, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        R = convert(Vector{TT}, [1, tableau.R∞])
        params = ParametersVPRKpSymmetric{DT,D}(equations, tableau, nullvec, Δt, NamedTuple{(:R,)}((R,)))

        # create cache dict
        caches = OldCacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*(S+1), params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpSymmetric(params, solver, iguess, caches)
    end

    function IntegratorVPRKpSymmetric(problem::Union{IODEProblem{DT},LODEProblem{DT}}, tableau, nullvec; kwargs...) where {DT}
        IntegratorVPRKpSymmetric{DT, ndims(problem)}(functions(problem), tableau, nullvec, timestep(problem); kwargs...)
    end
end


function Base.show(io::IO, int::IntegratorVPRKpSymmetric)
    print(io, "\nVariational Partitioned Runge-Kutta Integrator with Symmetric Projection and:\n")
    print(io, "   Timestep: $(int.timestep(problem))\n")
    print(io, "   Tableau:  $(description(int.tableau(method)))\n")
    print(io, "   $(string(int.tableau(method).q))")
    print(io, "   $(string(int.tableau(method).p))")
    # print(io, reference(int.tableau(method)))
end


function Integrators.get_internal_variables(int::IntegratorVPRKpSymmetric{DT,TT,D,S}) where {DT, TT, D, S}
    Q = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)
    λ = zeros(DT,D)

    solver = get_solver_status(int.solver)

    (Q=Q, P=P, V=V, F=F, λ=λ, solver=solver)
end


function initial_guess!(int::IntegratorVPRKpSymmetric{DT,TT}, sol::SolutionStepPODE{DT,TT},
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
                params::ParametersVPRKpSymmetric{DT,TT,D,S}) where {ST,DT,TT,D,S}

    local t₀::TT = solstep.t̄[1]
    local t₁::TT = solstep.t̄[1] + timestep(problem)
    local y1::ST
    local y2::ST

    # copy x to λ and q̄
    for k in 1:D
        λ[k] = x[D*S+k]
    end

    # compute U=λ
    U[1] .= λ
    U[2] .= λ

    # compute G=g(q,λ)
    for k in 1:D
        y1 = 0
        y2 = 0
        for j in 1:S
            y1 += tableau(method).q.b[j] * V[j][k]
            y2 += tableau(method).q.b̂[j] * V[j][k]
        end
        q[k] = solstep.q̄[1][k] + timestep(problem) * (y1 + y2) + timestep(problem) * (params.pparams[:R][1] * U[1][k] + params.pparams[:R][2] * U[2][k])
    end
    
    # compute G=g(q,λ)
    functions(problem).g(G[1], t₀, solstep.q̄[1], params.v̄, λ)
    functions(problem).g(G[2], t₁, q, v, λ)

    # compute p=ϑ(q)⋅p
    functions(problem).ϑ(p, t₁, q, v)
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSymmetric{DT,TT,D,S},
                caches::OldCacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.q̃, cache.p̃, cache.ṽ, cache.λ, cache.Q, cache.V, cache.U, cache.P, cache.F, cache.G, params)

    # compute b = - [P-AF-U]
    compute_rhs_vprk!(b, cache.P, cache.F, cache.G, params)

    # compute b = - [p-bF-G]
    compute_rhs_vprk_projection_p!(b, cache.p̃, cache.F, cache.G, D*(S+0), params)

    compute_rhs_correction!(b, cache.V, params)
end


function integrate_step!(int::IntegratorVPRKpSymmetric{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # reset solution
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

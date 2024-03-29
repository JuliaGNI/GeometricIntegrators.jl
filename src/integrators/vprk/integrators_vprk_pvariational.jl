
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpVariational = AbstractParametersVPRK{:vprk_pvariational}

function IntegratorCache(params::ParametersVPRKpVariational{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(D*S, true, 2*D; kwargs...)
end

function IntegratorCache{ST}(params::ParametersVPRKpVariational{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(D*S, true, 2*D; kwargs...)
end


"Variational Partitioned Runge-Kutta Integrator with Variational Projection."
struct IntegratorVPRKpVariational{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpVariational{DT,TT},
                ST  <: NonlinearSolver,
                PST <: NonlinearSolver,
                IT <: InitialGuessIODE{TT}} <: GeometricIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    caches::OldCacheDict{PPT}

    function IntegratorVPRKpVariational(params::ParametersVPRK{DT,TT,D,S},
                    pparams::ParametersVPRKpVariational{DT,TT,D,S},
                    solver::ST, projector::PST, iguess::IT, caches) where {DT,TT,D,S,ST,PST,IT}
        new{DT, TT, D, S, typeof(params), typeof(pparams), ST, PST, IT}(params, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKpVariational{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, nullvec, Δt::TT; R=[1,1]) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRK{DT,D}(equations, tableau, nullvec, Δt)

        # create projector params
        R  = convert(Vector{TT}, R)
        R1 = [R[1], zero(TT)]
        R2 = [zero(TT), R[2]]
        pparams = ParametersVPRKpVariational{DT,D}(equations, tableau, nullvec, Δt, NamedTuple{(:R,:R1,:R2)}((R,R1,R2)))

        # create cache dict
        caches = OldCacheDict(pparams)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create projector
        projector = create_nonlinear_solver(DT, 2*D, pparams, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKpVariational(params, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKpVariational(problem::Union{IODEProblem{DT},LODEProblem{DT}}, tableau, nullvec; kwargs...) where {DT}
        IntegratorVPRKpVariational{DT, ndims(problem)}(functions(problem), tableau, nullvec, timestep(problem); kwargs...)
    end
end


function Base.show(io::IO, int::IntegratorVPRKpVariational)
    print(io, "\nVariational Partitioned Runge-Kutta Integrator with Variational Projection and:\n")
    print(io, "   Timestep: $(int.timestep(problem))\n")
    print(io, "   Tableau:  $(description(int.tableau(method)))\n")
    print(io, "   $(string(int.tableau(method).q))")
    print(io, "   $(string(int.tableau(method).p))")
    # print(io, reference(int.tableau(method)))
end


function initialize!(int::IntegratorVPRKpVariational{DT}, sol::SolutionStepPODE{DT},
                                 cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    equation(int, :v̄)(sol.v, sol.t, sol.q)
    equation(int, :f̄)(sol.f, sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)

    # initialise projector
    equation(int, :g)(cache.G[1], sol.t, sol.q, sol.v, cache.λ)
    # cache.U[1] .= cache.λ
end


function initial_guess!(int::IntegratorVPRKpVariational{DT}, sol::SolutionStepPODE{DT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.history.q[2], sol.history.p[2], sol.history.v[2], sol.history.f[2],
                              sol.history.q[1], sol.history.p[1], sol.history.v[1], sol.history.f[1],
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            cache.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end

function initial_guess_projection!(int::IntegratorVPRKpVariational{DT}, sol::SolutionStepPODE{DT},
                                   cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    offset_q = 0
    offset_λ = ndims(int)

    for k in eachdim(int)
        cache.x̄[offset_q+k] = sol.q[k]
        cache.x̄[offset_λ+k] = 0
    end
end


function compute_projection!(
                x::Vector{ST}, q::AbstractVector{ST}, p::AbstractVector{ST}, v::AbstractVector{ST}, λ::AbstractVector{ST},
                U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpVariational{DT,TT,D,S}) where {ST,DT,TT,D,S}

    @assert D == length(q) == length(p) == length(v) == length(λ)
    @assert D == length(U[1]) == length(U[2])
    @assert D == length(G[1]) == length(G[2])

    # copy x to q, λ
    for k in 1:D
        q[k] = x[0*D+k]
        λ[k] = x[1*D+k]
    end

    # compute u=λ and g=∇ϑ(q)⋅λ
    U[1] .= 0
    U[2] .= λ
    # U[1] .= λ
    # U[2] .= 0

    functions(problem).g(G[1], solstep.t̄, q, v, λ)
    G[2] .= 0
    # G[1] .= 0
    # functions(problem).g(G[2], solstep.t̄, q, v, λ)

    # compute p=ϑ(q)
    functions(problem).ϑ(p, solstep.t̄, q, v)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
function residual!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpVariational{DT,TT,D,S},
                caches::OldCacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_projection!(x, cache.q̃, cache.p̃, cache.ṽ, cache.λ, cache.U, cache.G, params)

    # # compute b = - [q̄-q-U]
    for k in 1:D
        b[0*D+k] = - (cache.q̃[k] - solstep.q̄[k]) + timestep(problem) * params.pparams[:R][2] * cache.U[2][k]
        # b[0*D+k] = - (cache.q[k] - solstep.q̄[k])
    end

    # compute b = - [p̄-p-G]
    for k in 1:D
        b[1*D+k] = - (cache.p̃[k] - solstep.p̄[k])
        # b[1*D+k] = - (cache.p[k] - solstep.p̄[k]) + timestep(problem) * params.R[2] * cache.G[2][k]
    end
end


function integrate_step!(int::IntegratorVPRKpVariational{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int, sol, int.pparams.pparams[:R1], cache)

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

    # compute vector fields at internal stages
    components!(cache.x, cache.Q, cache.V, cache.P, cache.F, int.params)

    # compute unprojected solution
    update_solution!(int, sol, cache)

    # set time and solution for projection solver
    update_params!(int.pparams, sol)

    # set initial guess for projection
    initial_guess_projection!(int, sol, cache)

    # call projection solver
    solve!(cache.x̄, int.projector)

    # print solver status
    # print_solver_status(int.projector.status, int.projector.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.projector.status, int.projector.params)

    # compute projection vector fields
    compute_projection!(cache.x̄, cache.q̃, cache.p̃, cache.ṽ, cache.λ,
                                         cache.U, cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, sol, int.pparams.pparams[:R2], cache)
end

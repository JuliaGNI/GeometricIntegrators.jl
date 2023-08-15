
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKdegenerate = AbstractParametersVPRK{:vprk_degenerate}

function Integrators.IntegratorCache(params::ParametersVPRKdegenerate{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(D*S, true, D; kwargs...)
end

function Integrators.IntegratorCache{ST}(params::ParametersVPRKdegenerate{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(D*S, true, D; kwargs...)
end


"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKdegenerate{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKdegenerate{DT,TT},
                ST  <: NonlinearSolver,
                PST <: NonlinearSolver,
                IT  <: InitialGuessIODE{TT}} <: GeometricIntegratorVPRK{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    caches::OldCacheDict{PPT}

    function IntegratorVPRKdegenerate(params::ParametersVPRK{DT,TT,D,S},
                    pparams::ParametersVPRKdegenerate{DT,TT,D,S},
                    solver::ST, projector::PST, iguess::IT, caches) where {DT,TT,D,S,ST,PST,IT}
        new{DT, TT, D, S, typeof(params), typeof(pparams), ST, PST, IT}(params, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKdegenerate{DT,D}(equations::NamedTuple, tableau::PartitionedTableau{TT}, nullvec, Δt::TT; R=[1,1]) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # create solver params
        sparams = ParametersVPRK{DT,D}(equations, tableau, nullvec, Δt)

        # create projector params
        pparams = ParametersVPRKdegenerate{DT,D}(equations, tableau, nullvec, Δt, NamedTuple{(:V,:F)}((create_internal_stage_vector(DT,D,S), create_internal_stage_vector(DT,D,S))))

        # create cache dict
        caches = OldCacheDict(pparams)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, sparams, caches)

        # create projector
        projector = create_nonlinear_solver(DT, D, pparams, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKdegenerate(sparams, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKdegenerate(problem::Union{IODEProblem{DT},LODEProblem{DT}}, tableau, nullvec; kwargs...) where {DT}
        IntegratorVPRKdegenerate{DT, ndims(problem)}(functions(problem), tableau, nullvec, timestep(problem); kwargs...)
    end
end


function Base.show(io::IO, int::IntegratorVPRKdegenerate)
    print(io, "\nVariational Partitioned Runge-Kutta Integrator for Degenerate Lagrangians with:\n")
    print(io, "   Timestep: $(int.timestep(problem))\n")
    print(io, "   Tableau:  $(description(int.tableau(method)))\n")
    print(io, "   $(string(int.tableau(method).q))")
    print(io, "   $(string(int.tableau(method).p))")
    # print(io, reference(int.tableau(method)))
end


function initial_guess!(int::IntegratorVPRKdegenerate{DT}, sol::SolutionStepPODE{DT},
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

function initial_guess_projection!(int::IntegratorVPRKdegenerate{DT}, sol::SolutionStepPODE{DT},
                                   cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for k in eachdim(int)
        cache.x̄[k] = sol.q[k]
    end
end


function compute_solution!(x::Vector{ST}, q̄::Vector{ST}, v̄::Vector{ST}, p̄::Vector{ST},
                           params::ParametersVPRKdegenerate{DT,TT,D,S}) where {ST,DT,TT,D,S}

    @assert length(x) == length(q̄) == length(p̄)

    # copy x to q
    q̄ .= x

    # compute p̄ = ϑ(q)
    functions(problem).ϑ(p̄, solstep.t̄ + timestep(problem), q̄, v̄)
end


"Compute solution of degenerate symplectic partitioned Runge-Kutta methods."
function Integrators.residual!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,D,S},
                caches::OldCacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_solution!(x, cache.q̃, cache.ṽ, cache.p̃, params)

    # compute b = - [q̄-q-BV]
    for k in 1:div(D,2)
        b[k] = solstep.q̄[k] - cache.q̃[k]
        for i in 1:S
            b[k] += timestep(problem) * tableau(method).q.b[i] * params.pparams[:V][i][k]
        end
    end

    # compute b = - [p̄-p-BF]
    for k in 1:div(D,2)
        b[div(D,2)+k] = solstep.p̄[k] - cache.p̃[k]
        for i in 1:S
            b[div(D,2)+k] += timestep(problem) * tableau(method).p.b[i] * params.pparams[:F][i][k]
        end
    end
end


function integrate_step!(int::IntegratorVPRKdegenerate{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)
    update_params!(int.pparams, sol)

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

    # copy vector fields
    for i in eachstage(int)
        int.pparams.pparams[:V][i] .= cache.V[i]
        int.pparams.pparams[:F][i] .= cache.F[i]
    end

    # set initial guess for projection
    initial_guess_projection!(int, sol)

    # call solution solver
    solve!(cache.x̄, int.projector)

    # print solver status
    # print_solver_status(int.projector.status, int.projector.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.projector.status, int.projector.params)

    # compute projection vector fields
    compute_solution!(cache.x̄, sol.q, cache.ṽ, sol.p, int.pparams)
end

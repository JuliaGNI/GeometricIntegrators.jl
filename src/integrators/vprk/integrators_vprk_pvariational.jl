
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpVariational = AbstractParametersVPRK{:vprk_pvariational}


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpVariational{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpVariational{DT,TT},
                ST  <: NonlinearSolver{DT},
                PST <: NonlinearSolver{DT},
                IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    caches::CacheDict{PPT}

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
        caches = CacheDict(pparams)

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


function Integrators.initialize!(int::IntegratorVPRKpVariational{DT}, sol::SolutionStepPODE{DT},
                                 cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    sol.t̄ = sol.t - timestep(int)

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
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end

function initial_guess_projection!(int::IntegratorVPRKpVariational, sol::SolutionStepPODE)
    offset_q = 0
    offset_λ = ndims(int)

    for k in eachdim(int)
        int.projector.x[offset_q+k] = sol.q[k]
        int.projector.x[offset_λ+k] = 0
    end
end


function compute_projection!(
                x::Vector{ST}, q::SolutionVector{ST}, p::SolutionVector{ST}, v::SolutionVector{ST}, λ::SolutionVector{ST},
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

    params.equ.g(G[1], params.t̄, q, v, λ)
    G[2] .= 0
    # G[1] .= 0
    # params.equ.g(G[2], params.t̄, q, v, λ)

    # compute p=ϑ(q)
    params.equ.ϑ(p, params.t̄, q, v)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpVariational{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_projection!(x, cache.q̃, cache.p̃, cache.ṽ, cache.λ, cache.U, cache.G, params)

    # # compute b = - [q̄-q-U]
    for k in 1:D
        b[0*D+k] = - (cache.q̃[k] - params.q̄[k]) + params.Δt * params.pparams[:R][2] * cache.U[2][k]
        # b[0*D+k] = - (cache.q[k] - params.q̄[k])
    end

    # compute b = - [p̄-p-G]
    for k in 1:D
        b[1*D+k] = - (cache.p̃[k] - params.p̄[k])
        # b[1*D+k] = - (cache.p[k] - params.p̄[k]) + params.Δt * params.R[2] * cache.G[2][k]
    end
end


function Integrators.integrate_step!(int::IntegratorVPRKpVariational{DT,TT}, sol::SolutionStepPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int, sol, int.pparams.pparams[:R1], cache)

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

    # compute unprojected solution
    update_solution!(int, sol, cache)

    # set time and solution for projection solver
    update_params!(int.pparams, sol)

    # set initial guess for projection
    initial_guess_projection!(int, sol)

    # call projection solver
    solve!(int.projector)

    # print solver status
    print_solver_status(int.projector.status, int.projector.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.projector.status, int.projector.params)

    # compute projection vector fields
    compute_projection!(int.projector.x, cache.q̃, cache.p̃, cache.ṽ, cache.λ,
                                         cache.U, cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, sol, int.pparams.pparams[:R2], cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

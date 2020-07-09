
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKdegenerate = AbstractParametersVPRK{:vprk_degenerate}


"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKdegenerate{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKdegenerate{DT,TT},
                ST  <: NonlinearSolver{DT},
                PST <: NonlinearSolver{DT},
                IT  <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    caches::CacheDict{PPT}

    function IntegratorVPRKdegenerate(params::ParametersVPRK{DT,TT,D,S},
                    pparams::ParametersVPRKdegenerate{DT,TT,D,S},
                    solver::ST, projector::PST, iguess::IT, caches) where {DT,TT,D,S,ST,PST,IT}
        new{DT, TT, D, S, typeof(params), typeof(pparams), ST, PST, IT}(params, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKdegenerate{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT; R=[1,1]) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # create solver params
        sparams = ParametersVPRK{DT,D}(equations, tableau, Δt)

        # create projector params
        pparams = ParametersVPRKdegenerate{DT,D}(equations, tableau, Δt, NamedTuple{(:V,:F)}((create_internal_stage_vector(DT,D,S), create_internal_stage_vector(DT,D,S))))

        # create cache dict
        caches = CacheDict(pparams)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, sparams, caches)

        # create projector
        projector = create_nonlinear_solver(DT, D, pparams, caches)

        # create initial guess
        iguess = InitialGuessIODE{DT,D}(get_config(:ig_interpolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRKdegenerate(sparams, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKdegenerate(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRKdegenerate{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end



function initial_guess!(int::IntegratorVPRKdegenerate{DT}, sol::AtomicSolutionPODE{DT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])
        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end

function initial_guess_projection!(int::IntegratorVPRKdegenerate{DT}, sol::AtomicSolutionPODE{DT},
                                   cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for k in eachdim(int)
        int.projector.x[k] = sol.q[k]
    end
end


function compute_solution!(x::Vector{ST}, q̅::Vector{ST}, v̅::Vector{ST}, p̅::Vector{ST},
                           params::ParametersVPRKdegenerate{DT,TT,D,S}) where {ST,DT,TT,D,S}

    @assert length(x) == length(q̅) == length(p̅)

    # copy x to q
    q̅ .= x

    # compute p̅ = ϑ(q)
    params.equ[:ϑ](params.t̅ + params.Δt, q̅, v̅, p̅)
end


"Compute solution of degenerate symplectic partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_solution!(x, cache.q̃, cache.ṽ, cache.p̃, params)

    # compute b = - [q̅-q-BV]
    for k in 1:div(D,2)
        b[k] = params.q̅[k] - cache.q̃[k]
        for i in 1:S
            b[k] += params.Δt * params.tab.q.b[i] * params.pparams[:V][i][k]
        end
    end

    # compute b = - [p̅-p-BF]
    for k in 1:div(D,2)
        b[div(D,2)+k] = params.p̅[k] - cache.p̃[k]
        for i in 1:S
            b[div(D,2)+k] += params.Δt * params.tab.p.b[i] * params.pparams[:F][i][k]
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKdegenerate{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)
    update_params!(int.pparams, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset solution
    reset!(sol, timestep(int))

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

    # copy vector fields
    for i in eachstage(int)
        int.pparams.pparams[:V][i] .= cache.V[i]
        int.pparams.pparams[:F][i] .= cache.F[i]
    end

    # set initial guess for projection
    initial_guess_projection!(int, sol)

    # call solution solver
    solve!(int.projector)

    # print solver status
    print_solver_status(int.projector.status, int.projector.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.projector.status, int.projector.params)

    # compute projection vector fields
    compute_solution!(int.projector.x, sol.q, cache.ṽ, sol.p, int.pparams)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

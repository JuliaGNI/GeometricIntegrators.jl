
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpStandard = AbstractParametersVPRK{:vprk_pstandard}


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpStandard{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpStandard{DT,TT},
                ST  <: NonlinearSolver{DT},
                PST <: NonlinearSolver{DT},
                IT  <: InitialGuessIODE{DT,TT}} <: AbstractIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    caches::CacheDict{PPT}

    function IntegratorVPRKpStandard(params::ParametersVPRK{DT,TT,D,S},
                    pparams::ParametersVPRKpStandard{DT,TT,D,S},
                    solver::ST, projector::PST, iguess::IT, caches) where {DT,TT,D,S,ST,PST,IT}
        new{DT, TT, D, S, typeof(params), typeof(pparams), ST, PST, IT}(params, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKpStandard{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT,
                                        RU::Vector, RG::Vector; R∞::Int=tableau.R∞) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # create params
        params  = ParametersVPRK{DT,D}(equations, tableau, Δt)

        # create projector params
        RU  = Vector{TT}(RU)
        RG  = Vector{TT}(RG)
        RU1 = [RU[1], zero(TT)]
        RU2 = [zero(TT), R∞ * RU[2]]
        RG1 = [RG[1], zero(TT)]
        RG2 = [zero(TT), R∞ * RG[2]]

        pparams = ParametersVPRKpStandard{DT,D}(equations, tableau, Δt,
                    NamedTuple{(:RU, :RU1, :RU2, :RG, :RG1, :RG2)}((RU, RU1, RU2, RG, RG1, RG2)))

        # create cache dict
        caches = CacheDict(pparams)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create projector
        projector = create_nonlinear_solver(DT, 2*D, pparams, caches)

        # create initial guess
        iguess = InitialGuessIODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create integrator
        IntegratorVPRKpStandard(params, pparams, solver, projector, iguess, caches)
    end

    function IntegratorVPRKpStandard(equation::IODE{DT}, tableau, Δt, RU, RG; kwargs...) where {DT}
        IntegratorVPRKpStandard{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt, RU, RG; kwargs...)
    end
end

function IntegratorVPRKpStandard(args...)
    IntegratorVPRKpStandardConstructor(args..., [0,1]; R∞=1)
end

function IntegratorVPRKpSymplectic(args...)
    IntegratorVPRKpStandardConstructor(args..., [1,1])
end

function IntegratorVPRKpStandardConstructor(equation, tableau, Δt, R; R∞=tableau.R∞)
    IntegratorVPRKpStandard(equation, tableau, Δt, R, R; R∞=R∞)
end

function IntegratorVPRKpVariationalQ(args...)
    IntegratorVPRKpStandard(args..., [1,0], [0,1]; R∞=1)
end

function IntegratorVPRKpVariationalP(args...)
    IntegratorVPRKpStandard(args..., [0,1], [1,0]; R∞=1)
end


function Integrators.initialize!(int::IntegratorVPRKpStandard{DT}, sol::AtomicSolutionPODE{DT},
                                 cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v)(sol.t, sol.q, sol.v)
    equation(int, :f)(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)

    # initialise projector
    cache.U[1] .= cache.λ
    equation(int, :g)(sol.t, sol.q, cache.λ, cache.G[1])
end


function initial_guess!(int::IntegratorVPRKpStandard{DT}, sol::AtomicSolutionPODE{DT},
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

function initial_guess_projection!(int::IntegratorVPRKpStandard, sol::AtomicSolutionPODE)
    offset_q = 0
    offset_λ = ndims(int)

    for k in eachdim(int)
        int.projector.x[offset_q+k] = sol.q[k]
        int.projector.x[offset_λ+k] = 0
    end
end


function compute_projection!(
                x::Vector{ST}, q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpStandard{DT,TT,D,S}) where {ST,DT,TT,D,S}

    @assert D == length(q) == length(p) == length(λ)
    @assert D == length(U[1]) == length(U[2])
    @assert D == length(G[1]) == length(G[2])

    # copy x to q, λ
    for k in 1:D
        q[k] = x[0*D+k]
        λ[k] = x[1*D+k]
    end

    # compute u=λ and g=∇ϑ(q)⋅λ
    U[1] .= λ
    U[2] .= λ

    params.equ[:g](params.t̅, q, λ, G[1])
    G[2] .= G[1]

    # compute p̅=ϑ(q)
    params.equ[:ϑ](params.t̅, q, λ, p)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpStandard{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_projection!(x, cache.q̃, cache.p̃, cache.λ, cache.U, cache.G, params)

    # compute b = - [q̅-q-U]
    for k in 1:D
        b[0*D+k] = - (cache.q̃[k] - params.q̅[k]) + params.Δt * params.pparams[:RU][2] * cache.U[2][k]
    end

    # compute b = - [p̅-p-G]
    for k in 1:D
        b[1*D+k] = - (cache.p̃[k] - params.p̅[k]) + params.Δt * params.pparams[:RG][2] * cache.G[2][k]
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKpStandard{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int, sol, int.pparams.pparams[:RU1], int.pparams.pparams[:RG1], cache)

    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset cache
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
    compute_projection!(int.projector.x, cache.q̃, cache.p̃, cache.λ, cache.U, cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, sol, int.pparams.pparams[:RU2], int.pparams.pparams[:RG2], cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

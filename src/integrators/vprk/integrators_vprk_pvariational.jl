
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpVariational = AbstractParametersVPRK{:vprk_pvariational}


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpVariational{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpVariational{DT,TT},
                ST  <: NonlinearSolver{DT},
                PST <: NonlinearSolver{DT},
                IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}

    function IntegratorVPRKpVariational(params::ParametersVPRK{DT,TT,D,S},
                    pparams::ParametersVPRKpVariational{DT,TT,D,S},
                    solver::ST, projector::PST, iguess::IT, cache) where {DT,TT,D,S,ST,PST,IT}
        new{DT, TT, D, S, typeof(params), typeof(pparams), ST, PST, IT}(params, pparams, solver, projector, iguess, cache)
    end

    function IntegratorVPRKpVariational{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT; R=[1,1]) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRK{DT,D}(equations, tableau, Δt)

        # create projector params
        R  = convert(Vector{TT}, R)
        R1 = [R[1], zero(TT)]
        R2 = [zero(TT), R[2]]
        pparams = ParametersVPRKpVariational{DT,D}(equations, tableau, Δt, NamedTuple{(:R,:R1,:R2)}((R,R1,R2)))

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, params)

        # create projector
        projector = create_nonlinear_solver(DT, 2*D, pparams)

        # create initial guess
        iguess = InitialGuessPODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create cache
        cache = IntegratorCacheVPRK{DT,D,S}(true)

        # create integrator
        IntegratorVPRKpVariational(params, pparams, solver, projector, iguess, cache)
    end

    function IntegratorVPRKpVariational(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRKpVariational{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


function compute_projection!(
                x::Vector{ST}, q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpVariational{DT,TT,D,S}
            ) where {ST,DT,TT,D,S}

    @assert D == length(q) == length(p) == length(λ)
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

    params.equ[:g](params.t̅, q, λ, G[1])
    G[2] .= 0
    # G[1] .= 0
    # params.equ.g(params.t̅, q, λ, G[2])

    # compute p=ϑ(q)
    params.equ[:ϑ](params.t̅, q, λ, p)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
@generated function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpVariational{DT,TT,D,S}
            ) where {ST,DT,TT,D,S}

    cache = IntegratorCacheVPRK{ST, D, S}(true)

    function_stages = quote
        @assert length(x) == length(b)

        compute_projection!(x, $cache.q̃, $cache.p̃, $cache.λ, $cache.U, $cache.G, params)

        # # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q̃[k] - params.q̅[k]) + params.Δt * params.pparams[:R][2] * $cache.U[2][k]
            # b[0*D+k] = - ($cache.q[k] - params.q̅[k])
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p̃[k] - params.p̅[k])
            # b[1*D+k] = - ($cache.p[k] - params.p̅[k]) + params.Δt * params.R[2] * $cache.G[2][k]
        end
    end

    return function_stages
end


function initial_guess!(int::IntegratorVPRKpVariational, sol::AtomicSolutionPODE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = int.cache.ṽ[k]
        end
    end
end

function initial_guess_projection!(int::IntegratorVPRKpVariational, sol::AtomicSolutionPODE)
    offset_q = 0
    offset_λ = ndims(int)

    for k in eachdim(int)
        int.projector.x[offset_q+k] = sol.q[k]
        int.projector.x[offset_λ+k] = 0
    end
end

function Integrators.initialize!(int::IntegratorVPRKpVariational, sol::AtomicSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v)(sol.t, sol.q, sol.p, sol.v)
    equation(int, :f)(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)

    # initialise projector
    equation(int, :g)(sol.t, sol.q, int.cache.λ, int.cache.G[1])
    # cache.U[1] .= int.cache.λ
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKpVariational{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int, sol, int.pparams.pparams[:R1])

    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.params)

    # compute unprojected solution
    update_solution!(int, sol)

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
    compute_projection!(int.projector.x, int.cache.q̃, int.cache.p̃, int.cache.λ,
                                         int.cache.U, int.cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, sol, int.pparams.pparams[:R2])

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end


"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKdegenerate = AbstractParametersVPRK{:vprk_degenerate}


"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKdegenerate{DT, TT, D, S,
                PT  <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKdegenerate{DT,TT},
                ST  <: NonlinearSolver{DT},
                PST <: NonlinearSolver{DT},
                IT  <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}

    params::PT
    pparams::PPT
    solver::ST
    projector::PST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}

    function IntegratorVPRKdegenerate(params::ParametersVPRK{DT,TT,D,S},
                    pparams::ParametersVPRKdegenerate{DT,TT,D,S},
                    solver::ST, projector::PST, iguess::IT, cache) where {DT,TT,D,S,ST,PST,IT}
        new{DT, TT, D, S, typeof(params), typeof(pparams), ST, PST, IT}(params, pparams, solver, projector, iguess, cache)
    end

    function IntegratorVPRKdegenerate{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT; R=[1,1]) where {DT, TT, D}
        # get number of stages
        S = tableau.s

        # create solver params
        sparams = ParametersVPRK{DT,D}(equations, tableau, Δt)

        # create projector params
        pparams = ParametersVPRKdegenerate{DT,D}(equations, tableau, Δt, NamedTuple{(:V,:F)}((create_internal_stage_vector(DT,D,S), create_internal_stage_vector(DT,D,S))))

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, sparams)

        # create projector
        projector = create_nonlinear_solver(DT, D, pparams)

        # create initial guess
        iguess = InitialGuessPODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create cache
        cache = IntegratorCacheVPRK{DT,D,S}(true)

        # create integrator
        IntegratorVPRKdegenerate(sparams, pparams, solver, projector, iguess, cache)
    end

    function IntegratorVPRKdegenerate(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRKdegenerate{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end



function initial_guess!(int::IntegratorVPRKdegenerate, sol::AtomicSolutionPODE)
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

function initial_guess_projection!(int::IntegratorVPRKdegenerate, sol::AtomicSolutionPODE)
    for k in eachdim(int)
        int.projector.x[k] = sol.q[k]
    end
end


function compute_solution!(
                x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,D,S}
            ) where {ST,DT,TT,D,S}

    @assert length(x) == length(q̅) == length(p̅)

    # copy x to q
    q̅ .= x

    # compute p̅ = ϑ(q)
    params.equ[:ϑ](params.t̅ + params.Δt, q̅, p̅)

end


"Compute solution of degenerate symplectic partitioned Runge-Kutta methods."
@generated function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,D,S}
            ) where {ST,DT,TT,D,S}

    # create temporary vectors
    q̅ = zeros(ST,D)
    p̅ = zeros(ST,D)

    quote
        @assert length(x) == length(b)

        compute_solution!(x, $q̅, $p̅, params)

        # compute b = - [q̅-q-BV]
        for k in 1:div(D,2)
            b[k] = params.q̅[k] - $q̅[k]
            for i in 1:S
                b[k] += params.Δt * params.tab.q.b[i] * params.pparams[:V][i][k]
            end
        end

        # compute b = - [p̅-p-BF]
        for k in 1:div(D,2)
            b[div(D,2)+k] = params.p̅[k] - $p̅[k]
            for i in 1:S
                b[div(D,2)+k] += params.Δt * params.tab.p.b[i] * params.pparams[:F][i][k]
            end
        end
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKdegenerate{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)
    update_params!(int.pparams, sol)

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

    # copy vector fields
    for i in eachstage(int)
        int.pparams.pparams[:V][i] .= int.cache.V[i]
        int.pparams.pparams[:F][i] .= int.cache.F[i]
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
    compute_solution!(int.projector.x, sol.q, sol.p, int.pparams)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

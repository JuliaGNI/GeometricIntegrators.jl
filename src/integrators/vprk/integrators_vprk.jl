
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRK = AbstractParametersVPRK{:vprk}


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRK{DT, TT, D, S, PT <: ParametersVPRK{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}

    function IntegratorVPRK(params::ParametersVPRK{DT,TT,D,S}, solver::ST, iguess::IT, cache) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, cache)
    end

    function IntegratorVPRK{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRK{DT,D}(equations, tableau, Δt)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, params)

        # create initial guess
        iguess = InitialGuessPODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create cache
        cache = IntegratorCacheVPRK{DT,D,S}()

        # create integrator
        IntegratorVPRK(params, solver, iguess, cache)
    end

    function IntegratorVPRK{DT,D}(ϑ::Function, f::Function, g::Function, v::Function,
                    tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorVPRK{DT,D}(NamedTuple{(:ϑ,:f,:g,:v)}((ϑ, f, g, v)), tableau, Δt; kwargs...)
    end

    function IntegratorVPRK{DT,D}(ϑ::Function, f::Function, g::Function, v::Function, h::Function,
                    tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorVPRK{DT,D}(NamedTuple{(:ϑ,:f,:g,:v,:h)}((ϑ, f, g, v, h)), tableau, Δt; kwargs...)
    end

    function IntegratorVPRK(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRK{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


IntegratorVPRKpNone = IntegratorVPRK


function initial_guess!(int::IntegratorVPRK, sol::AtomicSolutionPODE)
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


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRK{DT,TT,D,S}
            ) where {ST,DT,TT,D,S}

    cache = IntegratorCacheVPRK{ST, D, S}()

    quote
        @assert length(x) == length(b)

        compute_stages!(x, $cache.Q, $cache.V, $cache.P, $cache.F, params)
        compute_rhs_vprk!(b, $cache.P, $cache.F, params)
        compute_rhs_vprk_correction!(b, $cache.V, params)
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
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

    # compute final update
    update_solution!(int, sol)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

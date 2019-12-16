
# using Printf

"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKdegenerate{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT
    cache::IntegratorCacheVPRK{DT,D,S}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersVPRKdegenerate(equ::ET, tab::TableauVPRK{TT}, Δt::TT,
                                  cache::IntegratorCacheVPRK{DT,D,S}) where {DT, TT, ET <: IODE{DT,TT}, D, S}

    @assert D == equ.d
    @assert S == tab.s

    q = zeros(DT,D)
    p = zeros(DT,D)

    ParametersVPRKdegenerate{DT, TT, ET, D, S}(equ, tab, Δt, cache, 0, q, p)
end


"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKdegenerate{DT, TT,
                SPT <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKdegenerate{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT  <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVPRK{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT

    cache::IntegratorCacheVPRK{DT,D,S}
end

function IntegratorVPRKdegenerate(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create solver params
    sparams = ParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = ParametersVPRKdegenerate(equation, tableau, Δt, cache)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheVPRK{DT,D,S}(true)

    # create integrator
    IntegratorVPRKdegenerate{DT, TT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess), D, S}(
                sparams, pparams, solver, projector, iguess, cache, q, p)
end


"Compute solution of degenerate symplectic partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    q̅ = zeros(ST,D)
    p̅ = zeros(ST,D)

    quote
        @assert length(x) == length(b)

        compute_solution!(x, $q̅, $p̅, params)

        # compute b = - [q̅-q-BV]
        for k in 1:div(D,2)
            b[k] = params.q[k] - $q̅[k]
            for i in 1:S
                b[k] += params.Δt * params.tab.q.b[i] * params.cache.V[k,i]
            end
        end

        # compute b = - [p̅-p-BF]
        for k in 1:div(D,2)
            b[div(D,2)+k] = params.p[k] - $p̅[k]
            for i in 1:S
                b[div(D,2)+k] += params.Δt * params.tab.p.b[i] * params.cache.F[k,i]
            end
        end
    end
end

function compute_solution!(
                x::Vector{ST}, q̅::Vector{ST}, p̅::Vector{ST},
                params::ParametersVPRKdegenerate{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    @assert length(x) == length(q̅) == length(p̅)

    # copy x to q
    q̅ .= x

    # compute p̅ = ϑ(q)
    params.equ.ϑ(params.t + params.Δt, q̅, p̅)

end

function initial_guess!(int::IntegratorVPRKdegenerate, sol::AtomisticSolutionPODE)
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

function initial_guess_projection!(int::IntegratorVPRKdegenerate, sol::AtomisticSolutionPODE)
    for k in eachdim(int)
        int.projector.x[k] = sol.q[k]
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKdegenerate{DT,TT}, sol::AtomisticSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.sparams, sol)
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
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.sparams)

    # compute unprojected solution
    update_solution!(int, sol)

    # set initial guess for projection
    initial_guess_projection!(int, sol)

    # call solution solver
    solve!(int.projector)

    # print solver status
    print_solver_status(int.projector.status, int.projector.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.projector.status, int.projector.params)

    # compute projection vector fields
    compute_projection!(int.projector.x, sol.q, sol.p, int.pparams)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

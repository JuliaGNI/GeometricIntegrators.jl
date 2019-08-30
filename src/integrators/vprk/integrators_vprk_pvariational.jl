
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpVariational{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    R::Vector{TT}
    R1::Vector{TT}
    R2::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
end

function ParametersVPRKpVariational(equ::ET, tab::TableauVPRK{TT}, Δt::TT, R::Vector) where {DT, TT, ET <: IODE{DT,TT}}
    R  = convert(Vector{TT}, R)
    R1 = [one(TT), zero(TT)]
    R2 = [zero(TT), one(TT)]

    q = zeros(DT, equ.d)
    p = zeros(DT, equ.d)

    ParametersVPRKpVariational{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, R, R1, R2, zero(TT), q, p)
end

function update_params!(params::ParametersVPRKpVariational, cache::AbstractIntegratorCacheVPRK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
end


function compute_projection!(
                x::Vector{ST}, q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpVariational{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

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

    params.equ.g(params.t, q, λ, G[1])
    G[2] .= 0
    # G[1] .= 0
    # params.equ.g(params.t, q, λ, G[2])

    # compute p=ϑ(q)
    params.equ.ϑ(params.t, q, λ, p)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpVariational{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
        @assert length(x) == length(b)

        compute_projection!(x, $cache.q, $cache.p, $cache.λ, $cache.U, $cache.G, params)

        # # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q[k] - params.q[k]) + params.Δt * params.R[2] * $cache.U[2][k]
            # b[0*D+k] = - ($cache.q[k] - params.q[k])
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p[k] - params.p[k])
            # b[1*D+k] = - ($cache.p[k] - params.p[k]) + params.Δt * params.R[2] * $cache.G[2][k]
        end
    end

    return function_stages
end

"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKpVariational{DT, TT,
                SPT <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpVariational{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRKwProjection{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT
end

function IntegratorVPRKpVariational(equation::ET, tableau::TableauVPRK{TT}, Δt::TT;
                                    R=[1,1]) where {DT, TT, ET <: IODE{DT,TT}}

    D = equation.d
    S = tableau.s

    # create solver params
    sparams = ParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = ParametersVPRKpVariational(equation, tableau, Δt, R)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, 2*D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRKpVariational{DT, TT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess)}(
            sparams, pparams, solver, projector, iguess)
end

equation(integrator::IntegratorVPRKpVariational) = integrator.sparams.equ
timestep(integrator::IntegratorVPRKpVariational) = integrator.sparams.Δt
tableau(integrator::IntegratorVPRKpVariational) = integrator.sparams.tab
nstages(integrator::IntegratorVPRKpVariational) = integrator.sparams.tab.s


function initialize!(int::IntegratorVPRKpVariational{DT,TT}, cache::IntegratorCacheVPRK) where {DT,TT}
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)

    # initialise projector
    equation(int).g(cache.t, cache.q, cache.λ, cache.G[1])
    # cache.U[1] .= cache.λ
end


function initial_guess!(int::IntegratorVPRKpVariational, cache::IntegratorCacheVPRK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in 1:ndims(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end

function initial_guess_projection!(int::IntegratorVPRKpVariational, cache::IntegratorCacheVPRK)
    offset_q = 0
    offset_λ = ndims(int)

    for k in 1:ndims(int)
        int.projector.x[offset_q+k] = cache.q[k]
        int.projector.x[offset_λ+k] = 0
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpVariational{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int ,cache, int.pparams.R1)

    # update nonlinear solver parameters from cache
    update_params!(int.sparams, cache)

    # compute initial guess
    initial_guess!(int, cache)

    # reset cache
    reset!(cache, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, cache.n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, cache.n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.P, cache.F, int.sparams)

    # compute unprojected solution
    update_solution!(int, cache)

    # set time and solution for projection solver
    update_params!(int.pparams, cache)

    # set initial guess for projection
    initial_guess_projection!(int, cache)

    # call projection solver
    solve!(int.projector)

    # print solver status
    print_solver_status(int.projector.status, int.projector.params, cache.n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.projector.status, int.projector.params, cache.n)

    # compute projection vector fields
    compute_projection!(int.projector.x, cache.q̃, cache.p̃, cache.λ, cache.U, cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, cache, int.pparams.R2)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end

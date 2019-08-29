
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpStandard{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    RU::Vector{TT}
    RG::Vector{TT}

    RU1::Vector{TT}
    RU2::Vector{TT}
    RG1::Vector{TT}
    RG2::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}

    function ParametersVPRKpStandard{DT,TT,ET,D,S}(equ, tab, Δt, RU, RG) where {DT,TT,ET,D,S}
        RU1 = [RU[1], zero(TT)]
        RU2 = [zero(TT), tab.R∞ * RU[2]]
        RG1 = [RG[1], zero(TT)]
        RG2 = [zero(TT), tab.R∞ * RG[2]]
        new(equ, tab, Δt, RU, RG, RU1, RU2, RG1, RG2, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end

function ParametersVPRKpStandard(equ::ET, tab::TableauVPRK{TT}, Δt::TT, R::Vector) where {DT, TT, ET <: IODE{DT,TT}}
    R = Vector{TT}(R)
    ParametersVPRKpStandard{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, R, R)
end

function ParametersVPRKpStandard(equ::ET, tab::TableauVPRK{TT}, Δt::TT, RU::Vector, RG::Vector) where {DT, TT, ET <: IODE{DT,TT}}
    RU = Vector{TT}(RU)
    RG = Vector{TT}(RG)
    ParametersVPRKpStandard{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, RU, RG)
end

function update_params!(params::ParametersVPRKpStandard, cache::AbstractIntegratorCacheVPRK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
end

function compute_projection!(
                x::Vector{ST}, q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpStandard{DT,TT,ET,D,S}
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
    U[1] .= λ
    U[2] .= λ

    params.equ.g(params.t, q, λ, G[1])
    G[2] .= G[1]

    # compute p̅=ϑ(q)
    params.equ.ϑ(params.t, q, λ, p)
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpStandard{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    quote
        @assert length(x) == length(b)

        compute_projection!(x, $cache.q, $cache.p, $cache.λ, $cache.U, $cache.G, params)

        # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q[k] - params.q[k]) + params.Δt * params.RU[2] * $cache.U[2][k]
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p[k] - params.p[k]) + params.Δt * params.RG[2] * $cache.G[2][k]
        end
    end
end

"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKpStandard{DT, TT,
                SPT <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpStandard{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT  <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRK{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT
end

function IntegratorVPRKpStandard(args...)
    IntegratorVPRKpStandardConstructor(args..., [0,1])
end

function IntegratorVPRKpSymplectic(args...)
    IntegratorVPRKpStandardConstructor(args..., [1,1])
end

function IntegratorVPRKpVariationalQ(args...)
    IntegratorVPRKpStandardConstructor(args..., [1,0], [0,1])
end

function IntegratorVPRKpVariationalP(args...)
    IntegratorVPRKpStandardConstructor(args..., [0,1], [1,0])
end

function IntegratorVPRKpStandardConstructor(equation, tableau, Δt, R)
    IntegratorVPRKpStandardConstructor(equation, tableau, Δt, R, R)
end

function IntegratorVPRKpStandardConstructor(equation::ET, tableau::TableauVPRK{TT}, Δt::TT,
                                    RU::Vector, RG::Vector) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create solver params
    sparams = ParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = ParametersVPRKpStandard(equation, tableau, Δt, RU, RG)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, 2*D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRKpStandard{DT, TT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess)}(
                sparams, pparams, solver, projector, iguess)
end

equation(integrator::IntegratorVPRKpStandard) = integrator.sparams.equ
timestep(integrator::IntegratorVPRKpStandard) = integrator.sparams.Δt
tableau(integrator::IntegratorVPRKpStandard) = integrator.sparams.tab
nstages(integrator::IntegratorVPRKpStandard) = integrator.sparams.tab.s


function create_integrator_cache(int::IntegratorVPRKpStandard{DT,TT}) where {DT,TT}
    IntegratorCacheVPRKwProjection(DT, TT, ndims(int), nstages(int))
end


function initialize!(int::IntegratorVPRKpStandard{DT,TT}, cache::IntegratorCacheVPRK) where {DT,TT}
    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)

    # initialise projector
    cache.U[1] .= cache.λ
    equation(int).g(cache.t, cache.q, cache.λ, cache.G[1])
end


function initial_guess!(int::IntegratorVPRKpStandard, cache::IntegratorCacheVPRK)
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

function initial_guess_projection!(int::IntegratorVPRKpStandard, cache::IntegratorCacheVPRK)
    offset_q = 0
    offset_λ = ndims(int)

    for k in 1:ndims(int)
        int.projector.x[offset_q+k] = cache.q[k]
        int.projector.x[offset_λ+k] = 0
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpStandard{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    update_solution!(cache.q, cache.U, int.pparams.RU1, int.pparams.Δt)
    update_solution!(cache.p, cache.G, int.pparams.RG1, int.pparams.Δt)

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
    update_solution!(cache.q, cache.V, int.sparams.tab.q.b, int.sparams.tab.q.b̂, int.sparams.Δt)
    update_solution!(cache.p, cache.F, int.sparams.tab.p.b, int.sparams.tab.p.b̂, int.sparams.Δt)

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
    update_solution!(cache.q, cache.U, int.pparams.RU2, int.pparams.Δt)
    update_solution!(cache.p, cache.G, int.pparams.RG2, int.pparams.Δt)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end

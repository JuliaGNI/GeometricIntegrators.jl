
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

    function ParametersVPRKpStandard{DT,TT,ET,D,S}(equ, tab, Δt, RU, RG, R∞) where {DT,TT,ET,D,S}
        RU1 = [RU[1], zero(TT)]
        RU2 = [zero(TT), R∞ * RU[2]]
        RG1 = [RG[1], zero(TT)]
        RG2 = [zero(TT), R∞ * RG[2]]
        new(equ, tab, Δt, RU, RG, RU1, RU2, RG1, RG2, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end

function ParametersVPRKpStandard(equ::ET, tab::TableauVPRK{TT}, Δt::TT, R::Vector, R∞=tab.R∞) where {DT, TT, ET <: IODE{DT,TT}}
    R = Vector{TT}(R)
    ParametersVPRKpStandard{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, R, R, R∞)
end

function ParametersVPRKpStandard(equ::ET, tab::TableauVPRK{TT}, Δt::TT, RU::Vector, RG::Vector, R∞=tab.R∞) where {DT, TT, ET <: IODE{DT,TT}}
    RU = Vector{TT}(RU)
    RG = Vector{TT}(RG)
    ParametersVPRKpStandard{DT, TT, ET, equ.d, tab.s}(equ, tab, Δt, RU, RG, R∞)
end

function update_params!(params::ParametersVPRKpStandard, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    params.t  = sol.t
    params.q .= sol.q
    params.p .= sol.p
end


"Variational partitioned Runge-Kutta integrator."
mutable struct IntegratorVPRKpStandard{DT, TT,
                SPT <: ParametersVPRK{DT,TT},
                PPT <: ParametersVPRKpStandard{DT,TT},
                SST <: NonlinearSolver{DT},
                STP <: NonlinearSolver{DT},
                IT  <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVPRKwProjection{DT,TT}

    sparams::SPT
    pparams::PPT

    solver::SST
    projector::STP
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}
end

function IntegratorVPRKpStandard(args...)
    IntegratorVPRKpStandardConstructor(args..., [0,1]; R∞=1)
end

function IntegratorVPRKpSymplectic(args...)
    IntegratorVPRKpStandardConstructor(args..., [1,1])
end

function IntegratorVPRKpVariationalQ(args...)
    IntegratorVPRKpStandardConstructor(args..., [1,0], [0,1]; R∞=1)
end

function IntegratorVPRKpVariationalP(args...)
    IntegratorVPRKpStandardConstructor(args..., [0,1], [1,0]; R∞=1)
end

function IntegratorVPRKpStandardConstructor(equation, tableau, Δt, R; R∞=tableau.R∞)
    IntegratorVPRKpStandardConstructor(equation, tableau, Δt, R, R; R∞=R∞)
end

function IntegratorVPRKpStandardConstructor(equation::ET, tableau::TableauVPRK{TT}, Δt::TT,
                                    RU::Vector, RG::Vector; R∞::Int=tableau.R∞) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create solver params
    sparams = ParametersVPRK(equation, tableau, Δt)

    # create projector params
    pparams = ParametersVPRKpStandard(equation, tableau, Δt, RU, RG, R∞)

    # create nonlinear solver
    solver = create_nonlinear_solver(DT, D*S, sparams)

    # create projector
    projector = create_nonlinear_solver(DT, 2*D, pparams)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheVPRK{DT,D,S}(true)

    # create integrator
    IntegratorVPRKpStandard{DT, TT, typeof(sparams), typeof(pparams), typeof(solver), typeof(projector), typeof(iguess), D, S}(
                sparams, pparams, solver, projector, iguess, cache)
end

equation(integrator::IntegratorVPRKpStandard) = integrator.sparams.equ
timestep(integrator::IntegratorVPRKpStandard) = integrator.sparams.Δt
tableau(integrator::IntegratorVPRKpStandard) = integrator.sparams.tab


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

    cache = IntegratorCacheVPRK{ST, D, S}(true)

    quote
        @assert length(x) == length(b)

        compute_projection!(x, $cache.q̃, $cache.p̃, $cache.λ, $cache.U, $cache.G, params)

        # compute b = - [q̅-q-U]
        for k in 1:D
            b[0*D+k] = - ($cache.q̃[k] - params.q[k]) + params.Δt * params.RU[2] * $cache.U[2][k]
        end

        # compute b = - [p̅-p-G]
        for k in 1:D
            b[1*D+k] = - ($cache.p̃[k] - params.p[k]) + params.Δt * params.RG[2] * $cache.G[2][k]
        end
    end
end


function initialize!(int::IntegratorVPRKpStandard, sol::AtomicSolutionPODE)
    sol.t̅ = sol.t - timestep(int)

    equation(int).v(sol.t, sol.q, sol.p, sol.v)
    equation(int).f(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)

    # initialise projector
    int.cache.U[1] .= int.cache.λ
    equation(int).g(sol.t, sol.q, int.cache.λ, int.cache.G[1])
end


function initial_guess!(int::IntegratorVPRKpStandard, sol::AtomicSolutionPODE)
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

function initial_guess_projection!(int::IntegratorVPRKpStandard, sol::AtomicSolutionPODE)
    offset_q = 0
    offset_λ = ndims(int)

    for k in eachdim(int)
        int.projector.x[offset_q+k] = sol.q[k]
        int.projector.x[offset_λ+k] = 0
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpStandard{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # add perturbation for next time step to solution
    # (same vector field as previous time step)
    project_solution!(int, sol, int.pparams.RU1, int.pparams.RG1)

    # update nonlinear solver parameters from cache
    update_params!(int.sparams, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
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
    compute_projection!(int.projector.x, int.cache.q̃, int.cache.p̃, int.cache.λ, int.cache.U, int.cache.G, int.pparams)

    # add projection to solution
    project_solution!(int, sol, int.pparams.RU2, int.pparams.RG2)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

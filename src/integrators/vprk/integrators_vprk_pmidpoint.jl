
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpMidpoint{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    R::Vector{TT}

    t̅::TT
    q̅::Vector{DT}
    p̅::Vector{DT}
end

function ParametersVPRKpMidpoint(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    R = convert(Vector{TT}, [1, tab.R∞])
    q̅ = zeros(DT, ndims(equ))
    p̅ = zeros(DT, ndims(equ))

    ParametersVPRKpMidpoint{DT, TT, ET, ndims(equ), tab.s}(equ, tab, Δt, R, zero(TT), q̅, p̅)
end


@generated function compute_projection_vprk!(x::Vector{ST},
                q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpMidpoint{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    # create temporary vectors
    q̃  = zeros(ST,D)
    # qm = zeros(ST,D)

    quote
        local t₀::TT = params.t̅
        local t₁::TT = params.t̅ + params.Δt
        local tₘ::TT = (t₀+t₁)/2
        local y::ST

        # copy x to λ and q̅
        for k in 1:D
            q[k] = x[D*(S+0)+k]
            λ[k] = x[D*(S+1)+k]
        end

        # compute U=λ
        U[1] .= λ
        U[2] .= λ

        # compute G=g(q,λ)
        for k in 1:D
            y = 0
            for j in 1:S
                y += params.tab.q.b[j] * V[j][k]
            end
            $q̃[k] = params.q̅[k] + 0.5 * params.Δt * y + params.Δt * params.R[1] * U[1][k]
            # $qm[k] = params.q̅[k] + 0.5 * params.Δt * y + 0.5 * params.Δt * params.R[1] * U[k,1] + 0.5 * params.Δt * params.R[2] * U[k,2]
        end

        # println("q̃mid = ", $q̃)
        # println("qmid = ", $qm)

        params.equ.g(tₘ, $q̃, λ, G[1])
        G[2] .= G[1]

        # compute p=ϑ(q)
        params.equ.ϑ(t₁, q, λ, p)
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpMidpoint{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    fcache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    function_stages = quote
        compute_stages!(x, $pcache.q, $pcache.p, $pcache.λ, $fcache.Q, $fcache.V, $pcache.U, $fcache.P, $fcache.F, $pcache.G, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $fcache.P, $fcache.F, $pcache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_projection_p!(b, $pcache.p, $fcache.F, $pcache.G, D*(S+0), params)

        # compute b = - [q-bV-U]
        compute_rhs_vprk_projection_q!(b, $pcache.q, $fcache.V, $pcache.U, D*(S+1), params)

        compute_rhs_vprk_correction!(b, $fcache.V, params)
    end

    return function_stages
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpMidpoint{DT, TT, PT <: ParametersVPRKpMidpoint{DT,TT},
                                       ST <: NonlinearSolver{DT},
                                       IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRKwProjection{DT,TT}

    params::PT
    solver::ST
    iguess::IT
end

function IntegratorVPRKpMidpoint(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpMidpoint(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*(S+2), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRKpMidpoint{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(params, solver, iguess)
end

equation(integrator::IntegratorVPRKpMidpoint) = integrator.params.equ
timestep(integrator::IntegratorVPRKpMidpoint) = integrator.params.Δt
tableau(integrator::IntegratorVPRKpMidpoint) = integrator.params.tab
nstages(integrator::IntegratorVPRKpMidpoint) = integrator.params.tab.s


function initialize!(int::IntegratorVPRKpMidpoint{DT,TT}, cache::IntegratorCacheVPRK) where {DT,TT}
    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end


function initial_guess!(int::IntegratorVPRKpMidpoint{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in 1:ndims(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end

    evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                          cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                          cache.q̃, one(TT))
    for k in 1:ndims(int)
        int.solver.x[ndims(int)*(nstages(int)+0)+k] = cache.q̃[k]
    end
    for k in 1:ndims(int)
        int.solver.x[ndims(int)*(nstages(int)+1)+k] = 0
    end
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpMidpoint{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, cache)

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

    # compute vector fields at internal stages and projection vector fields
    compute_stages!(int.solver.x,
                    cache.q̃, cache.p̃, cache.λ,
                    cache.Q, cache.V, cache.U,
                    cache.P, cache.F, cache.G, int.params)

    # compute unprojected solution
    update_solution!(int, cache)

    # add projection to solution
    project_solution!(int, cache, int.params.R)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end

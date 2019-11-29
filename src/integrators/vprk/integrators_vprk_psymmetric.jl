
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpSymmetric{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    R::Vector{TT}

    t̅::TT
    q̅::Vector{DT}
    p̅::Vector{DT}
end

function ParametersVPRKpSymmetric(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    R = convert(Vector{TT}, [1, tab.R∞])
    q̅ = zeros(DT, ndims(equ))
    p̅ = zeros(DT, ndims(equ))

    ParametersVPRKpSymmetric{DT, TT, ET, ndims(equ), tab.s}(equ, tab, Δt, R, zero(TT), q̅, p̅)
end


function compute_projection_vprk!(x::Vector{ST},
                q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpSymmetric{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    local t₀::TT = params.t̅
    local t₁::TT = params.t̅ + params.Δt

    # copy x to λ and q̅
    for k in 1:D
        q[k] = x[D*(S+0)+k]
        λ[k] = x[D*(S+1)+k]
    end

    # compute U=λ
    U[1] .= λ
    U[2] .= λ

    # compute G=g(q,λ)
    params.equ.g(t₀, params.q̅, λ, G[1])
    params.equ.g(t₁, q, λ, G[2])

    # compute p=ϑ(q)⋅p
    params.equ.ϑ(t₁, q, λ, p)
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSymmetric{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    scache = NonlinearFunctionCacheVPRK{ST}(D,S)
    pcache = NonlinearFunctionCacheVPRKprojection{ST}(D,S)

    quote
        compute_stages!(x, $pcache.q, $pcache.p, $pcache.λ, $scache.Q, $scache.V, $pcache.U, $scache.P, $scache.F, $pcache.G, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $scache.P, $scache.F, $pcache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_projection_p!(b, $pcache.p, $scache.F, $pcache.G, D*(S+0), params)

        # compute b = - [q-bV-U]
        compute_rhs_vprk_projection_q!(b, $pcache.q, $scache.V, $pcache.U, D*(S+1), params)

        compute_rhs_vprk_correction!(b, $scache.V, params)
    end
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpSymmetric{DT, TT, PT <: ParametersVPRKpSymmetric{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}} <: AbstractIntegratorVPRKwProjection{DT,TT}

    params::PT
    solver::ST
    iguess::IT
end

function IntegratorVPRKpSymmetric(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpSymmetric(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*(S+2), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorVPRKpSymmetric{DT, TT, typeof(params), typeof(solver), typeof(iguess)}(params, solver, iguess)
end

equation(integrator::IntegratorVPRKpSymmetric) = integrator.params.equ
timestep(integrator::IntegratorVPRKpSymmetric) = integrator.params.Δt
tableau(integrator::IntegratorVPRKpSymmetric) = integrator.params.tab
nstages(integrator::IntegratorVPRKpSymmetric) = integrator.params.tab.s


function initial_guess!(int::IntegratorVPRKpSymmetric{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
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
function integrate_step!(int::IntegratorVPRKpSymmetric{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, cache)

    # compute initial guess
    initial_guess!(int, cache)

    # reset cache
    reset!(cache, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

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

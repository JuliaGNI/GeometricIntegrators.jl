
"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
mutable struct AbstractParametersVSPARK{IT,DT,TT,D,S,R,FT,PT,UT,GT,ϕT} <: Parameters{DT,TT}
    f_f::FT
    f_p::PT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    t_q::CoefficientsARK{TT}
    t_p::CoefficientsARK{TT}
    t_q̃::CoefficientsPRK{TT}
    t_p̃::CoefficientsPRK{TT}
    t_λ::CoefficientsMRK{TT}
    t_ω::Matrix{TT}
    d_v::Vector{TT}

    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    function AbstractParametersVSPARK{IT,DT,TT,D,S,R,FT,PT,UT,GT,ϕT}(f_f, f_p, f_u, f_g, f_ϕ, Δt, t_q, t_p, t_q̃, t_p̃, t_λ, t_ω, d_v) where {IT,DT,TT,D,S,R,FT,PT,UT,GT,ϕT}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new(f_f, f_p, f_u, f_g, f_ϕ, Δt,
            t_q, t_p, t_q̃, t_p̃, t_λ, t_ω, d_v,
            zero(TT), q, p, λ)
    end
end


function update_params!(params::AbstractParametersVSPARK, cache::IntegratorCacheVSPARK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
    params.λ .= cache.λ
end


equation(int::AbstractIntegratorVSPARK) = int.equation
timestep(int::AbstractIntegratorVSPARK) = int.params.Δt
tableau(int::AbstractIntegratorVSPARK) = int.tableau
nstages(int::AbstractIntegratorVSPARK) = int.tableau.s
pstages(int::AbstractIntegratorVSPARK) = int.tableau.r


function create_integrator_cache(int::AbstractIntegratorVSPARK{DT,TT}) where {DT,TT}
    IntegratorCacheVSPARK{DT, TT, ndims(int), nstages(int), pstages(int)}()
end


function initialize!(int::AbstractIntegratorVSPARK, cache::IntegratorCacheVSPARK)
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end


function initial_guess!(int::AbstractIntegratorVSPARK, cache::IntegratorCacheVSPARK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+3] = cache.ṽ[k]
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if int.params.t_λ.c[1] == 0
        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = cache.λ[k]
        end
    end

    if isdefined(tableau(int), :d)
        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end


"Integrate an implicit DAE with a partitioned additive Runge-Kutta integrator for variational systems."
function integrate_step!(int::AbstractIntegratorVSPARK{DT,TT}, cache::IntegratorCacheVSPARK{DT,TT}) where {DT,TT}
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

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # compute final update
    update_solution!(cache.q, cache.Vi, int.params.t_q.b, timestep(int))
    update_solution!(cache.p, cache.Fi, int.params.t_p.b, timestep(int))

    # compute projection
    update_solution!(cache.q, cache.Up, int.params.t_q.β, timestep(int))
    update_solution!(cache.p, cache.Gp, int.params.t_p.β, timestep(int))
    # TODO # update_multiplier!(cache.λ, cache.Λp, int.params.t_λ.b)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end

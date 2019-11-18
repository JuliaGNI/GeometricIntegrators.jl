
"Integrate an implicit DAE with a specialised partitioned additive Runge-Kutta integrator."
function integrate_step!(int::AbstractIntegratorSPARK{DT,TT}, cache::IntegratorCacheSPARK{DT,TT}) where {DT,TT}
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
    update_solution!(cache.q, cache.qₑᵣᵣ, cache.Vi, int.params.tab.q.b, timestep(int))
    update_solution!(cache.p, cache.pₑᵣᵣ, cache.Fi, int.params.tab.p.b, timestep(int))

    # compute projection
    update_solution!(cache.q, cache.qₑᵣᵣ, cache.Up, int.params.tab.q.β, timestep(int))
    update_solution!(cache.p, cache.pₑᵣᵣ, cache.Gp, int.params.tab.p.β, timestep(int))
    # TODO # update_multiplier!(cache.λ, cache.Λp, int.params.tab.λ.b)

    # copy solution to initial guess
    update!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f)

    # take care of periodic solutions
    cut_periodic_solution!(cache, equation(int).periodicity)
end

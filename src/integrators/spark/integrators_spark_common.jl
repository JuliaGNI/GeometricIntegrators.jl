
"Integrate an implicit DAE with a specialised partitioned additive Runge-Kutta integrator."
function integrate_step!(int::AbstractIntegratorSPARK{DT,TT}, sol::AtomisticSolutionPDAE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

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
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.Vi, int.params.t_q.b, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Fi, int.params.t_p.b, timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, int.cache.Up, int.params.t_q.β, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Gp, int.params.t_p.β, timestep(int))
    # TODO # update_multiplier!(sol.λ, int.cache.Λp, int.params.t_λ.b)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

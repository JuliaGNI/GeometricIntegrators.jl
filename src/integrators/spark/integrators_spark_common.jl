
function Integrators.get_internal_variables(int::AbstractIntegratorSPARK{DT,TT,D,S,R}) where {DT, TT, D, S, R}
    Qi = create_internal_stage_vector(DT, D, S)
    Pi = create_internal_stage_vector(DT, D, S)
    Vi = create_internal_stage_vector(DT, D, S)
    Φi = create_internal_stage_vector(DT, D, S)

    Qp = create_internal_stage_vector(DT, D, R)
    Pp = create_internal_stage_vector(DT, D, R)
    Λp = create_internal_stage_vector(DT, D, R) 
    Φp = create_internal_stage_vector(DT, D, R)

    solver = get_solver_status(int.solver)

    (Qi=Qi, Pi=Pi, Vi=Vi, Φi=Φi, Qp=Qp, Pp=Pp, Λp=Λp, Φp=Φp, solver=solver)
end

function update_solution!(int::AbstractIntegratorSPARK{DT,TT}, sol::SolutionStepPDAE{DT,TT},
                          cache::IntegratorCacheSPARK{DT}=int.caches[DT]) where {DT,TT}
    # compute final update
    update_solution!(sol.q, sol.q̃, cache.Vi, int.params.tab.q.b, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.Fi, int.params.tab.p.b, timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, cache.Up, int.params.tab.q.β, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.Gp, int.params.tab.p.β, timestep(int))
    update_multiplier!(sol.λ, cache.Λp, int.params.tab.λ.b)
end


function Integrators.integrate_step!(int::AbstractIntegratorSPARK{DT,TT}, sol::SolutionStepPDAE{DT,TT},
                                     cache::IntegratorCacheSPARK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset cache
    reset!(sol)

    # call nonlinear solver
    solve!(int.solver)

    # check_jacobian(int.solver)
    # print_jacobian(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache, int.params)

    # compute final update
    update_solution!(int, sol, cache)

    # copy internal stage variables
    sol.internal.Qi .= cache.Qi
    sol.internal.Pi .= cache.Pi
    sol.internal.Vi .= cache.Vi
    sol.internal.Φi .= cache.Φi

    sol.internal.Qp .= cache.Qp
    sol.internal.Pp .= cache.Pp
    sol.internal.Λp .= cache.Λp
    sol.internal.Φp .= cache.Φp

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end

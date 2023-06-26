
function Integrators.get_internal_variables(method::AbstractSPARKMethod, problem::AbstractSPARKProblem{DT,TT}) where {DT,TT}
    S = nstages(method)
    R = pstages(method)
    D = ndims(problem)
    
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


function update_solution!(
    solstep::SolutionStepPDAE{DT,TT}, 
    problem::AbstractSPARKProblem{DT,TT},
    method::AbstractSPARKMethod, 
    caches::CacheDict) where {DT,TT}

    # compute final update
    # update_solution!(solstep.q, solstep.q̃, caches[DT].Vi, tableau(method).q.b, timestep(problem))
    # update_solution!(solstep.p, solstep.p̃, caches[DT].Fi, tableau(method).p.b, timestep(problem))
    update!(solstep.q, solstep.q̃, caches[DT].Vi, tableau(method).q.b, tableau(method).q.b̂, timestep(problem))
    update!(solstep.p, solstep.p̃, caches[DT].Fi, tableau(method).p.b, tableau(method).p.b̂, timestep(problem))

    # compute projection
    # update_solution!(solstep.q, solstep.q̃, caches[DT].Up, tableau(method).q.β, timestep(problem))
    # update_solution!(solstep.p, solstep.p̃, caches[DT].Gp, tableau(method).p.β, timestep(problem))
    update!(solstep.q, solstep.q̃, caches[DT].Up, tableau(method).q.β, tableau(method).q.β̂, timestep(problem))
    update!(solstep.p, solstep.p̃, caches[DT].Gp, tableau(method).p.β, tableau(method).p.β̂, timestep(problem))
    update_multiplier!(solstep.λ, caches[DT].Λp, tableau(method).λ.b)
end


function Integrators.integrate_step!(
    solstep::SolutionStepPDAE{DT},
    problem::AbstractSPARKProblem{DT},
    method::AbstractSPARKMethod,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT}

    # call nonlinear solver
    solve!(caches[DT].x, (b,x) -> residual!(b, x, solstep, problem, method, caches), solver)

    # check_jacobian(int.solver)
    # print_jacobian(int.solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute vector fields at internal stages
    components!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    update_solution!(solstep, problem, method, caches)

    # copy internal stage variables
    # TODO: reactivate
    # solstep.internal.Qi .= caches[DT].Qi
    # solstep.internal.Pi .= caches[DT].Pi
    # solstep.internal.Vi .= caches[DT].Vi
    # solstep.internal.Φi .= caches[DT].Φi

    # solstep.internal.Qp .= caches[DT].Qp
    # solstep.internal.Pp .= caches[DT].Pp
    # solstep.internal.Λp .= caches[DT].Λp
    # solstep.internal.Φp .= caches[DT].Φp
end

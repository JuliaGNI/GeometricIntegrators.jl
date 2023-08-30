
function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::Union{IODEProblem,LODEProblem},
    method::Union{HPImidpoint,HPItrapezoidal}, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    # get cache and dimension
    cache = caches[DT]
    D = ndims(problem)
    A = nparams(method)

    # compute initial guess for solution q(n+1)
    initialguess!(solstep.t, cache.q, cache.p, solstep, problem, iguess)

    # copy initial guess to solution vector
    cache.x[1:D] .= cache.q
    cache.x[D+1:D+A] .= 0
end

function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem{DT,TT},LODEProblem{DT,TT}},
    method::Union{HPImidpoint,HPItrapezoidal}, 
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, (b,x) -> residual!(b, x, solstep, problem, method, caches), solver)

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    components!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    solstep.q .= caches[DT].q
    solstep.p .= caches[DT].p
end

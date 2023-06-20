
function Integrators.initial_guess!(
    solstep::SolutionStepPDAE{DT}, 
    problem::Union{HDAEProblem,PDAEProblem},
    method::Union{HSPARKMethod,PSPARKMethod}, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    cache = caches[DT]

    for i in 1:nstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q.c[i], cache.Qi[i], cache.Pi[i], cache.Vi[i], cache.Fi[i], solstep, problem, iguess; nowarn = true)

        for k in 1:ndims(problem)
            cache.x[2*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qi[i][k] - solstep.q̄[k]) / timestep(problem)
            cache.x[2*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pi[i][k] - solstep.p̄[k]) / timestep(problem)
        end
    end

    for i in 1:pstages(method)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q̃.c[i], cache.Qp[i], cache.Pp[i], cache.Vp[i], cache.Fp[i], solstep, problem, iguess; nowarn = true)

        for k in 1:ndims(problem)
            cache.x[2*ndims(problem)*nstages(method)+3*(ndims(problem)*(i-1)+k-1)+1] = (cache.Qp[i][k] - solstep.q̄[k]) / timestep(problem)
            cache.x[2*ndims(problem)*nstages(method)+3*(ndims(problem)*(i-1)+k-1)+2] = (cache.Pp[i][k] - solstep.p̄[k]) / timestep(problem)
            cache.x[2*ndims(problem)*nstages(method)+3*(ndims(problem)*(i-1)+k-1)+3] = 0
        end
    end

    # TODO: Check indices !!!
    # if isdefined(tableau(method), :λ) && tableau(method).λ.c[1] == 0
    #     for k in 1:ndims(problem)
    #         cache.x[2*ndims(problem)*nstages(method)+3*(k-1)+3] = solstep.λ[k]
    #     end
    # end
end

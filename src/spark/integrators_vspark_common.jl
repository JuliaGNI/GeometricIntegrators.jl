
function initial_guess!(sol, history, params, int::GeometricIntegrator{<:VPARK,<:Union{IDAEProblem,LDAEProblem}})
    # get cache for internal stages
    local C = cache(int)

    for i in 1:nstages(int)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        soltmp = (
            t = history.t[1] + timestep(int) * tableau(int).q.c[i],
            q = cache(int).Qi[i],
            p = cache(int).Pi[i],
            v = cache(int).Vi[i],
            f = cache(int).Fi[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[3*(ndims(int)*(i-1)+k-1)+1] = (C.Qi[i][k] - sol.q[k]) / timestep(int)
            C.x[3*(ndims(int)*(i-1)+k-1)+2] = (C.Pi[i][k] - sol.p[k]) / timestep(int)
            C.x[3*(ndims(int)*(i-1)+k-1)+3] =  C.Vi[i][k]
        end
    end

    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        soltmp = (
            t = history.t[1] + timestep(int) * tableau(int).q̃.c[i],
            q = cache(int).Qp[i],
            p = cache(int).Pp[i],
            v = cache(int).Vp[i],
            f = cache(int).Fp[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:ndims(int)
            C.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (C.Qp[i][k] - sol.q[k]) / timestep(int)
            C.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (C.Pp[i][k] - sol.p[k]) / timestep(int)
            C.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    # TODO: Check indices !!!
    # if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
    #     for k in 1:ndims(int)
    #         C.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = sol.λ[k]
    #     end
    # end

    if hasnullvector(method(int))
        for k in 1:ndims(int)
            C.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(method(int))+k] = 0
        end
    end
end


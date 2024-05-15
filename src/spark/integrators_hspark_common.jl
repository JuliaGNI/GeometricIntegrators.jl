function initial_guess!(int::GeometricIntegrator{<:Union{HSPARKMethod,PSPARKMethod},<:Union{HDAEProblem,PDAEProblem}})
    # get cache for internal stages
    local C = cache(int)
    local sol = current(solstep(int))

    for i in 1:nstages(int)
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(sol.t + timestep(int) * (tableau(int).q.c[i] - 1), C.Qi[i], C.Pi[i], C.Vi[i], C.Fi[i], solstep(int), problem(int), iguess(int); nowarn = true)

        for k in 1:ndims(int)
            C.x[2*(ndims(int)*(i-1)+k-1)+1] = (C.Qi[i][k] - sol.q[k]) / timestep(int)
            C.x[2*(ndims(int)*(i-1)+k-1)+2] = (C.Pi[i][k] - sol.p[k]) / timestep(int)
        end
    end

    for i in 1:pstages(method(int))
        # TODO: initialguess! should take two timesteps for c[i] of q and p tableau
        initialguess!(sol.t + timestep(int) * (tableau(int).q̃.c[i] - 1), C.Qp[i], C.Pp[i], C.Vp[i], C.Fp[i], solstep(int), problem(int), iguess(int); nowarn = true)

        for k in 1:ndims(int)
            C.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (C.Qp[i][k] - sol.q[k]) / timestep(int)
            C.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (C.Pp[i][k] - sol.p[k]) / timestep(int)
            C.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    # TODO: Check indices !!!
    # if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
    #     for k in 1:ndims(int)
    #         C.x[2*ndims(int)*nstages(int)+3*(k-1)+3] = solstep(int).λ[k]
    #     end
    # end
end

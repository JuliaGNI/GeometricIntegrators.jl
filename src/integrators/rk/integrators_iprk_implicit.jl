
solversize(problem::AbstractProblemIODE, method::IPRK) = ndims(problem) * nstages(method)

function Base.show(io::IO, int::GeometricIntegrator{<:IPRK,<:AbstractProblemIODE})
    print(io, "\nPartitioned Runge-Kutta Integrator for Implicit Equations with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end

function Cache{ST}(problem::AbstractProblemIODE, method::IPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IPRKCache{ST,D,S,solversize(problem, method)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::IPRK) = IPRKCache{ST,ndims(problem),nstages(tableau(method)),solversize(problem, method)}


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:IPRK,<:AbstractProblemIODE})
    # get cache for nonlinear solution vector
    local x = nlsolution(int)

    # compute initial guess for internal stages
    for i in eachstage(int)
        soltmp = (
            t=history.t[1] + timestep(int) * tableau(int).p.c[i],
            q=cache(int).Q[i],
            p=cache(int).P[i],
            v=cache(int).V[i],
            f=cache(int).F[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(int)
        offset = ndims(int) * (i - 1)
        for k in 1:ndims(int)
            x[offset+k] = cache(int).P[i][k] - sol.p[k]
            for j in eachstage(int)
                x[offset+k] -= timestep(int) * tableau(int).p.a[i, j] * cache(int).F[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IPRK,<:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to V
    for i in eachindex(C.V)
        for k in eachindex(C.V[i])
            C.V[i][k] = x[ndims(int)*(i-1)+k]
        end
    end

    # compute Q = q + Δt A V
    for i in eachindex(C.Q)
        for k in eachindex(C.Q[i])
            y1 = y2 = zero(ST)
            for j in eachindex(C.V)
                y1 += tableau(int).q.a[i, j] * C.V[j][k]
                y2 += tableau(int).q.â[i, j] * C.V[j][k]
            end
            C.Q[i][k] = sol.q[k] + timestep(int) * (y1 + y2)
        end
    end

    # compute ϑ(Q,V) and f(Q,V)
    for i in eachindex(C.P, C.F)
        equations(int).ϑ(C.P[i], sol.t + timestep(int) * (tableau(int).q.c[i] - 1), C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], sol.t + timestep(int) * (tableau(int).p.c[i] - 1), C.Q[i], C.V[i], params)
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IPRK,<:AbstractProblemIODE}) where {ST}
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # get cache for internal stages
    local C = cache(int, ST)

    # compute b = - [(Y-AV), (Z-AF)]
    for i in eachstage(int)
        for k in eachindex(C.P[i])
            z1 = z2 = zero(ST)
            for j in eachstage(int)
                z1 += tableau(int).p.a[i, j] * C.F[j][k]
                z2 += tableau(int).p.â[i, j] * C.F[j][k]
            end
            b[ndims(int)*(i-1)+k] = C.P[i][k] - sol.p[k] - timestep(int) * (z1 + z2)
        end
    end
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:IPRK,<:AbstractProblemIODE}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol.q, cache(int, DT).V, tableau(int).q, timestep(int))
    update!(sol.p, cache(int, DT).F, tableau(int).p, timestep(int))
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:IPRK,<:AbstractProblemIODE})
    # call nonlinear solver
    solve!(solver(int), nlsolution(int), (sol, params, int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end

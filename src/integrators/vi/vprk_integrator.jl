
description(::GeometricIntegrator{<:VPRKMethod}) = "Variational Partitioned Runge-Kutta Integrator"

solversize(problem::AbstractProblemIODE, method::VPRKMethod) =
    ndims(problem) * nstages(method)


function internal_variables(method::VPRKMethod, problem::AbstractProblemIODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = ndims(problem)

    Q = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)
    Y = create_internal_stage_vector(DT, D, S)
    Z = create_internal_stage_vector(DT, D, S)

    # solver = get_solver_status(int.solver)

    (Q=Q, P=P, V=V, F=F, Y=Y, Z=Z)#, solver=solver)
end


function copy_internal_variables!(solstep::SolutionStep, cache::VPRKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :P) && copyto!(internal(solstep).P, cache.P)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :F) && copyto!(internal(solstep).F, cache.F)
    haskey(internal(solstep), :Y) && copyto!(internal(solstep).Y, cache.Y)
    haskey(internal(solstep), :Z) && copyto!(internal(solstep).Z, cache.Z)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:VPRK})
    for i in eachstage(int)
        soltmp = (
            t=history.t[1] + timestep(int) * tableau(int).q.c[i],
            q=cache(int).Q[i],
            p=cache(int).P[i],
            v=cache(int).V[i],
            f=cache(int).F[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
        for k in eachindex(cache(int).V[i])
            cache(int).x[ndims(int)*(i-1)+k] = cache(int).V[i][k]
        end
        # println("  t = $(solstep.t̄ + timestep(problem) * tableau(method).q.c[i]),",
        #         "  q̄ = $(solstep.q̄), v̄ = $(solstep.v̄), ",
        #         "  q = $(cache.Q[i]), v = $(cache.V[i])")
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRK}) where {ST}
    # get cache for internal stages
    local Q = cache(int, ST).Q
    local P = cache(int, ST).P
    local V = cache(int, ST).V
    local F = cache(int, ST).F

    # copy x to V
    for i in eachindex(V)
        for k in eachindex(V[i])
            V[i][k] = x[ndims(int)*(i-1)+k]
        end
    end

    # compute Q = q + Δt A V
    for i in eachindex(Q)
        @assert ndims(int) == length(Q[i]) == length(V[i])
        for k in eachindex(Q[i])
            y1 = y2 = zero(ST)
            for j in eachindex(V)
                y1 += tableau(int).q.a[i, j] * V[j][k]
                y2 += tableau(int).q.â[i, j] * V[j][k]
            end
            Q[i][k] = sol.q[k] + timestep(int) * (y1 + y2)
        end
    end

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in eachindex(P, F)
        tᵢ = sol.t + timestep(int) * (tableau(int).p.c[i] - 1)
        equations(int).ϑ(P[i], tᵢ, Q[i], V[i], params)
        equations(int).f(F[i], tᵢ, Q[i], V[i], params)
    end
end


function residual_solution!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRK}) where {ST}
    # get cache for previous solution and internal stages
    local P = cache(int, ST).P
    local F = cache(int, ST).F

    # compute b = - [(P-p-AF)]
    for i in eachindex(P)
        for k in eachindex(P[i])
            z1 = z2 = zero(ST)
            for j in eachindex(F)
                z1 += tableau(int).p.a[i, j] * F[j][k]
                z2 += tableau(int).p.â[i, j] * F[j][k]
            end
            b[ndims(int)*(i-1)+k] = -(P[i][k] - sol.p[k]) + timestep(int) * (z1 + z2)
        end
    end
end


function residual_correction!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRK}) where {ST}
    # get cache for internal stages
    local V = cache(int, ST).V
    local μ = cache(int, ST).μ

    local sl::Int = div(nstages(int) + 1, 2)

    if hasnullvector(int)
        # compute μ
        for k in eachindex(μ)
            μ[k] = tableau(int).p.b[sl] / nullvector(int)[sl] * b[ndims(int)*(sl-1)+k]
        end

        # replace equation for Pₗ with constraint on V
        for k in eachindex(μ)
            b[ndims(int)*(sl-1)+k] = 0
            for i in eachindex(nullvector(int))
                b[ndims(int)*(sl-1)+k] += V[i][k] * nullvector(int)[i]
            end
        end

        # modify P₁, ..., Pₛ except for Pₗ
        for i in eachindex(nullvector(int))
            if i ≠ sl
                z = nullvector(int)[i] / tableau(int).p.b[i]
                for k in eachindex(μ)
                    b[ndims(int)*(i-1)+k] -= z * μ[k]
                end
            end
        end
    end
end


# Compute stages of variational partitioned Runge-Kutta methods.
function residual!(b::AbstractVector, sol, params, int::GeometricIntegrator{<:VPRK})
    residual_solution!(b, sol, params, int)
    residual_correction!(b, sol, params, int)
end


# Compute stages of Variational Partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:VPRK}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, int::GeometricIntegrator{<:VPRK}, DT)
    # compute final update
    update!(sol.q, cache(int, DT).V, tableau(int).q, timestep(int))
    update!(sol.p, cache(int, DT).F, tableau(int).p, timestep(int))
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:VPRK}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:VPRK,<:AbstractProblemIODE})
    # call nonlinear solver
    solve!(solver(int), nlsolution(int), (sol, params, int))

    # check_jacobian(solver(int))
    # print_jacobian(solver(int))

    # print solver status
    # println(status(solver(int)))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver(int))))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end


solversize(problem::AbstractProblemIODE, method::IPRK) = ndims(problem) * nstages(method)

function Base.show(io::IO, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE})
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
    IPRKCache{ST, D, S, solversize(problem, method)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::IPRK) = IPRKCache{ST, ndims(problem), nstages(tableau(method)), solversize(problem, method)}


function initial_guess!(int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE})
    # get cache for nonlinear solution vector and internal stages
    local x = nlsolution(int)
    local Q = cache(int).Q
    local P = cache(int).P
    local V = cache(int).V
    local F = cache(int).F

    # compute initial guess for internal stages
    for i in eachstage(int)
        initialguess!(solstep(int).t̄ + timestep(int) * tableau(int).p.c[i], Q[i], P[i], V[i], F[i], solstep(int), problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(int)
        offset = ndims(int)*(i-1)
        for k in 1:ndims(int)
            x[offset+k] = P[i][k] - cache(int).p̄[k]
            for j in eachstage(int)
                x[offset+k] -= timestep(int) * tableau(int).p.a[i,j] * F[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local q̄ = cache(int, ST).q̄
    local Q = cache(int, ST).Q
    local P = cache(int, ST).P
    local V = cache(int, ST).V
    local F = cache(int, ST).F

    # copy x to V
    for i in eachindex(V)
        for k in eachindex(V[i])
            V[i][k] = x[ndims(int)*(i-1) + k]
        end
    end

    # compute Q = q + Δt A V
    for i in eachindex(Q)
        for k in eachindex(Q[i])
            y1 = y2 = zero(ST)
            for j in eachindex(V)
                y1 += tableau(int).q.a[i,j] * V[j][k]
                y2 += tableau(int).q.â[i,j] * V[j][k]
            end
            Q[i][k] = q̄[k] + timestep(int) * (y1 + y2)
        end
    end

    # compute ϑ(Q,V) and f(Q,V)
    for i in eachindex(P,F)
        equations(int).ϑ(P[i], solstep(int).t̄ + timestep(int) * tableau(int).q.c[i], Q[i], V[i], parameters(solstep(int)))
        equations(int).f(F[i], solstep(int).t̄ + timestep(int) * tableau(int).p.c[i], Q[i], V[i], parameters(solstep(int)))
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local p̄ = cache(int, ST).p̄
    local P = cache(int, ST).P
    local F = cache(int, ST).F

    # compute b = - [(Y-AV), (Z-AF)]
    for i in eachstage(int)
        for k in eachindex(P[i])
            z1 = z2 = zero(ST)
            for j in eachstage(int)
                z1 += tableau(int).p.a[i,j] * F[j][k]
                z2 += tableau(int).p.â[i,j] * F[j][k]
            end
            b[ndims(int)*(i-1)+k] = P[i][k] - p̄[k] - timestep(int) * (z1 + z2)
        end
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    update!(solstep(int), cache(int, DT).V, cache(int, DT).F, tableau(int), timestep(int))
end


function integrate_step!(int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(nlsolution(int), int)
end

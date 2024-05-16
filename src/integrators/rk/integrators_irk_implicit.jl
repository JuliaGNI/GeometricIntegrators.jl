@doc raw"""
Implicit Runge-Kutta integrator cache.

### Fields

* `q̄`: solution at previous timestep
* `p̄`: momentum at previous timestep
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Θ`: internal stages of one-form ``\vartheta``
* `F`: internal stages of force field
"""
struct IRKimplicitCache{DT,D,S} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}

    q̄::Vector{DT}
    p̄::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Θ::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function IRKimplicitCache{DT,D,S}(method::IRK) where {DT,D,S}
        if implicit_update(method)
            x = zeros(DT, D*(S+1))
        else
            x = zeros(DT, D*S)
        end

        q̄ = zeros(DT,D)
        p̄ = zeros(DT,D)

        q = zeros(DT,D)
        v = zeros(DT,D)
        θ = zeros(DT,D)
        f = zeros(DT,D)

        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Θ = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)
        
        new(x, q̄, p̄, q, v, θ, f, Q, V, Θ, F)
    end
end

nlsolution(cache::IRKimplicitCache) = cache.x

function reset!(cache::IRKimplicitCache, t, q, p, λ = missing)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

function Cache{ST}(problem::AbstractProblemIODE, method::IRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IRKimplicitCache{ST,D,S}(method; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::IRK) = IRKimplicitCache{ST, ndims(problem), nstages(tableau(method))}


function solversize(problem::AbstractProblemIODE, method::IRK)
    n = ndims(problem) * nstages(method)

    if implicit_update(method)
        n += ndims(problem)
    end

    return n
end


function internal_variables(method::IRK, problem::AbstractProblemIODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = ndims(problem)

    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    Θ = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)

    # solver = get_solver_status(int.solver)

    (Q=Q, V=V, Θ=Θ, F=F)#, solver=solver)
end

function copy_internal_variables(solstep::SolutionStep, cache::IRKimplicitCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :Θ) && copyto!(internal(solstep).Θ, cache.Θ)
    haskey(internal(solstep), :F) && copyto!(internal(solstep).F, cache.F)
end


function Base.show(io::IO, int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE})
    print(io, "\nRunge-Kutta Integrator for Implicit Equations with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


function initial_guess!(int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE})
    # get cache for nonlinear solution vector and internal stages
    local x = nlsolution(int)
    local Q = cache(int).Q
    local Θ = cache(int).Θ
    local V = cache(int).V
    local F = cache(int).F

    # compute initial guess for internal stages
    for i in eachstage(int)
        initialguess!(solstep(int).t̄ + timestep(int) * tableau(int).c[i], Q[i], Θ[i], V[i], F[i], solstep(int), problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(int)
        offset = ndims(int)*(i-1)
        for k in 1:ndims(int)
            x[offset+k] = Θ[i][k] - solstep(int).p̄[k]
            for j in eachstage(int)
                x[offset+k] -= timestep(int) * tableau(int).a[i,j] * F[j][k]
            end
        end
    end

    # compute initial guess for solution
    if implicit_update(int)
        initialguess!(solstep(int).t, cache(int).q, cache(int).θ, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

        offset = ndims(int) * nstages(int)

        for k in 1:ndims(int)
            x[offset + k] = cache(int).q[k]
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to V
    for i in eachindex(C.V)
        for k in eachindex(C.V[i])
            C.V[i][k] = x[ndims(int)*(i-1) + k]
        end
    end

    # copy x to q
    if implicit_update(int)
        for k in eachindex(C.q)
            C.q[k] = x[ndims(int) * nstages(int) + k]
        end
    end

    # compute Q = q + Δt A V
    for i in eachindex(C.Q)
        for k in eachindex(C.Q[i])
            y1 = y2 = zero(ST)
            for j in eachindex(C.V)
                y1 += tableau(int).a[i,j] * C.V[j][k]
                y2 += tableau(int).â[i,j] * C.V[j][k]
            end
            C.Q[i][k] = sol.q[k] + timestep(int) * (y1 + y2)
        end
    end

    # compute Θ = ϑ(Q) and F = f(Q,V)
    for i in eachindex(C.Θ, C.F)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        equations(int).ϑ(C.Θ[i], tᵢ, C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], tᵢ, C.Q[i], C.V[i], params)
    end

    # compute θ = ϑ(q)
    if implicit_update(int)
        equations(int).ϑ(C.θ, sol.t, C.q, C.v, params)
    end
end


# Compute stages of implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, sol, int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # compute b for internal stages
    for i in eachindex(C.Θ)
        for k in eachindex(C.Θ[i])
            y1 = y2 = zero(ST)
            for j in eachindex(C.F)
                y1 += tableau(int).a[i,j] * C.F[j][k]
                y2 += tableau(int).â[i,j] * C.F[j][k]
            end
            b[ndims(int)*(i-1) + k] = C.Θ[i][k] - sol.p[k] - timestep(int) * (y1 + y2)
        end
    end

    # compute b for update
    if implicit_update(int)
        for k in eachindex(C.θ)
            y1 = y2 = zero(ST)
            for j in eachindex(C.F)
                y1 += tableau(int).b[j] * C.F[j][k]
                y2 += tableau(int).b̂[j] * C.F[j][k]
            end
            b[ndims(int) * nstages(int) + k] = C.θ[k] - sol.p[k] - timestep(int) * (y1 + y2)
        end
    end
end


# Compute stages of implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), sol...)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), sol...)

    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol.q, cache(int, DT).V, tableau(int), timestep(int))
    update!(sol.p, cache(int, DT).F, tableau(int), timestep(int))
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:IRK, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end

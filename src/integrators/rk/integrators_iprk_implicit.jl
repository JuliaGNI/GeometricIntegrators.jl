
solversize(problem::AbstractProblemIODE, method::IPRK) = ndims(problem) * nstages(method)


function Base.show(io::IO, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE})
    print(io, "\nPartitioned Runge-Kutta Integrator for Implicit Equations with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end


@doc raw"""
Implicit partitioned Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of q
* `p̃`: initial guess of p
* `ṽ`: initial guess of v
* `f̃`: initial guess of f
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of q
* `P`: internal stages of p
* `V`: internal stages of v
* `F`: internal stages of f
* `Y`: vector field of internal stages of q
* `Z`: vector field of internal stages of p
"""
struct IPRKimplicitCache{ST,D,S} <: PODEIntegratorCache{ST,D}
    x::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    function IPRKimplicitCache{ST,D,S}() where {ST,D,S}
        # create solver vector
        x = zeros(ST, D*S)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)

        new(x, Q, P, V, F)
    end
end

nlsolution(cache::IPRKimplicitCache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::IPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IPRKimplicitCache{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::IPRK) = IPRKimplicitCache{ST, ndims(problem), nstages(tableau(method))}


function initial_guess!(int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE})
    # get cache for internal stages
    local x = cache(int).x
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
            x[offset+k] = P[i][k] - solstep(int).p̄[k]
            for j in eachstage(int)
                x[offset+k] -= timestep(int) * tableau(int).p.a[i,j] * F[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local Q = cache(int, ST).Q
    local P = cache(int, ST).P
    local V = cache(int, ST).V
    local F = cache(int, ST).F
    local D = ndims(int)

    for i in eachstage(int)
        for k in 1:D
            # copy y to V
            V[i][k] = x[D*(i-1)+k]

            # compute Q
            y1 = y2 = zero(ST)
            for j in eachstage(int)
                y1 += tableau(int).q.a[i,j] * V[j][k]
                y2 += tableau(int).q.â[i,j] * V[j][k]
            end
            Q[i][k] = solstep(int).q̄[k] + timestep(int) * (y1 + y2)
        end

        # compute ϑ(Q,V) and f(Q,V)
        equations(int).ϑ(P[i], solstep(int).t̄ + timestep(int) * tableau(int).q.c[i], Q[i], V[i])
        equations(int).f(F[i], solstep(int).t̄ + timestep(int) * tableau(int).p.c[i], Q[i], V[i])
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE}) where {ST}
    # get cache for internal stages
    local P = cache(int, ST).P
    local F = cache(int, ST).F
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute b = - [(Y-AV), (Z-AF)]
    for i in eachstage(int)
        for k in 1:D
            z1 = z2 = zero(ST)
            for j in eachstage(int)
                z1 += tableau(int).p.a[i,j] * F[j][k]
                z2 += tableau(int).p.â[i,j] * F[j][k]
            end
            b[D*(i-1)+k] = P[i][k] - solstep(int).p̄[k] - timestep(int) * (z1 + z2)
        end
    end
end


function integrate_step!(int::GeometricIntegrator{<:IPRK, <:AbstractProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute vector fields at internal stages
    components!(nlsolution(int), int)

    # compute final update
    update!(solstep(int), cache(int).V, cache(int).F, tableau(int), timestep(int))
end

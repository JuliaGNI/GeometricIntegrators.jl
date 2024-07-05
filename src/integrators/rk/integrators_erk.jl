@doc raw"""
Explicit Runge-Kutta method solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
"""
abstract type ERKMethod <: RKMethod end

isexplicit(method::Union{ERKMethod, Type{<:ERKMethod}}) = true
isimplicit(method::Union{ERKMethod, Type{<:ERKMethod}}) = false


"""
Explicit Runge-Kutta Method

```
ERK(tableau)
```
"""
struct ERK{TT <: Tableau} <: ERKMethod
    tableau::TT
end

initmethod(method::ERKMethod) = ERK(method)


function Base.show(io::IO, int::GeometricIntegrator{<:ERK})
    print(io, "\nExplicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(RungeKutta.description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(tableau(int)))
end


"Explicit Runge-Kutta integrator cache."
struct ERKCache{DT,D,S} <: ODEIntegratorCache{DT,D}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}

    function ERKCache{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        new(Q, V)
    end
end

function Cache{ST}(problem::AbstractProblemODE, method::ERK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    ERKCache{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemODE, method::ERK) = ERKCache{ST, ndims(problem), nstages(tableau(method))}

function internal_variables(method::ERK, problem::AbstractProblemODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = ndims(problem)

    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)

    (Q=Q, V=V)
end

function copy_internal_variables!(solstep::SolutionStep, cache::ERKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
end


function components!(_, sol, params, int::GeometricIntegrator{<:ERK, <:AbstractProblemODE})
    # obtain cache
    local Q = cache(int).Q
    local V = cache(int).V

    # compute internal stages
    for i in eachstage(int)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        for k in eachindex(Q[i], V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += tableau(int).a[i,j] * V[j][k]
            end
            Q[i][k] = sol.q[k] + timestep(int) * yᵢ
        end
        equations(int).v(V[i], tᵢ, Q[i], params)
    end
end

function update!(sol, params, _, int::GeometricIntegrator{<:ERK})
    # compute vector field at internal stages
    components!(nothing, sol, params, int)

    # compute final update
    update!(sol.q, cache(int).V, tableau(int), timestep(int))
end

function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ERK, <:AbstractProblemODE})
    # compute final update
    update!(sol, params, nothing, int)
end

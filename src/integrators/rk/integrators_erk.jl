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


function integrate_step!(int::GeometricIntegrator{<:ERK, <:AbstractProblemODE})
    # obtain cache
    local Q = cache(int).Q
    local V = cache(int).V

    # compute internal stages
    for i in eachstage(int)
        tᵢ = solstep(int).t̄ + timestep(int) * tableau(int).c[i]
        for k in eachindex(Q[i], V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += tableau(int).a[i,j] * V[j][k]
            end
            Q[i][k] = solstep(int).q̄[k] + timestep(int) * yᵢ
        end
        equations(int).v(V[i], tᵢ, Q[i], parameters(solstep(int)))
    end

    # compute final update
    update!(solstep(int), V, tableau(int), timestep(int))
end

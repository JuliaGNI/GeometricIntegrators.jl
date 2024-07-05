@doc raw"""
Explicit partitioned Runge-Kutta method solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (t_i, Q_{n,i}, P_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} .
\end{aligned}
```
Usually we are interested in Hamiltonian systems, where
```math
\begin{aligned}
V_{n,i} &= \dfrac{\partial H}{\partial p} (Q_{n,i}, P_{n,i}) , &
F_{n,i} &= - \dfrac{\partial H}{\partial q} (Q_{n,i}, P_{n,i}) , 
\end{aligned}
```
and tableaus satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```
"""
abstract type EPRKMethod <: PRKMethod end

isexplicit(method::Union{EPRKMethod, Type{<:EPRKMethod}}) = true
isimplicit(method::Union{EPRKMethod, Type{<:EPRKMethod}}) = false


"""
Explicit Partitioned Runge-Kutta Method

```
EPRK(tableau)
```
"""
struct EPRK{TT <: PartitionedTableau} <: EPRKMethod
    tableau::TT
end

initmethod(method::EPRKMethod) = EPRK(method)

function Base.show(io::IO, int::GeometricIntegrator{<:EPRK})
    print(io, "\nExplicit Partitioned Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end


"Explicit Runge-Kutta integrator cache."
struct EPRKCache{DT,D,S} <: PODEIntegratorCache{DT,D}
    Q::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}
    P::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}

    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    function EPRKCache{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector_with_zero(DT, D, S)
        P = create_internal_stage_vector_with_zero(DT, D, S)

        Y = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)

        new(Q, P, Y, Z, V, F)
    end
end

function Cache{ST}(problem::EquationProblem, method::EPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    EPRKCache{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::EPRK) = EPRKCache{ST, ndims(problem), nstages(tableau(method))}


function internal_variables(method::EPRK, problem::AbstractProblemPODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = ndims(problem)

    Q = create_internal_stage_vector_with_zero(DT, D, S)
    P = create_internal_stage_vector_with_zero(DT, D, S)

    Y = create_internal_stage_vector(DT, D, S)
    Z = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)

    (Q=Q, P=P, V=V, F=F, Y=Y, Z=Z)
end

function copy_internal_variables!(solstep::SolutionStep, cache::EPRKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :P) && copyto!(internal(solstep).P, cache.P)
    haskey(internal(solstep), :Y) && copyto!(internal(solstep).Y, cache.Y)
    haskey(internal(solstep), :Z) && copyto!(internal(solstep).Z, cache.Z)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :F) && copyto!(internal(solstep).F, cache.F)
end


# Compute Q stages of explicit partitioned Runge-Kutta methods.
function compute_stage_q!(sol, params, int::GeometricIntegrator{<:EPRK, <:AbstractProblemPODE}, i, jmax, t)
    # obtain cache
    local C = cache(int)

    fill!(C.Y[i], 0)
    for k in eachindex(C.Y[i])
        for j in 1:jmax
            C.Y[i][k] += timestep(int) * tableau(int).q.a[i,j] * C.V[j][k]
        end
    end
    for k in eachindex(C.Y[i])
        C.Q[i][k] = sol.q[k] + C.Y[i][k]
    end
    equations(int).f(C.F[i], t, C.Q[i], C.P[jmax], params)
end

# Compute P stages of explicit partitioned Runge-Kutta methods.
function compute_stage_p!(sol, params, int::GeometricIntegrator{<:EPRK, <:AbstractProblemPODE}, i, jmax, t)
    # obtain cache
    local C = cache(int)

    fill!(C.Z[i], 0)
    for k in eachindex(C.Z[i])
        for j in 1:jmax
            C.Z[i][k] += timestep(int) * tableau(int).p.a[i,j] * C.F[j][k]
        end
    end
    for k in eachindex(C.Z[i])
        C.P[i][k] = sol.p[k] + C.Z[i][k]
    end
    equations(int).v(C.V[i], t, C.Q[jmax], C.P[i], params)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:EPRK, <:AbstractProblemPODE})
    # store previous solution
    cache(int).Q[0] .= sol.q
    cache(int).P[0] .= sol.p

    # compute internal stages
    for i in eachstage(int)
        tqᵢ = sol.t + timestep(int) * (tableau(int).q.c[i] - 1)
        tpᵢ = sol.t + timestep(int) * (tableau(int).p.c[i] - 1)

        if tableau(int).q.a[i,i] ≠ 0 && tableau(int).p.a[i,i] ≠ 0
            error("This is an implicit method!")
        elseif tableau(int).q.a[i,i] ≠ 0
            compute_stage_p!(sol, params, int, i, i-1, tpᵢ)
            compute_stage_q!(sol, params, int, i, i, tqᵢ)
        elseif tableau(int).p.a[i,i] ≠ 0
            compute_stage_q!(sol, params, int, i, i-1, tqᵢ)
            compute_stage_p!(sol, params, int, i, i, tpᵢ)
        else
            compute_stage_q!(sol, params, int, i, i-1, tqᵢ)
            compute_stage_p!(sol, params, int, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update!(sol.q, cache(int).V, tableau(int).q, timestep(int))
    update!(sol.p, cache(int).F, tableau(int).p, timestep(int))
end

@doc raw"""
Implicit partitioned Runge-Kutta method solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (t_i, Q_{n,i}, P_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} ,
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

For implicit equations like [`IODE`](@ref)s and [`LODE`](@ref)s this method is solving the system
```math
\begin{aligned}
P_{n,i} &= \vartheta (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} ,
\end{aligned}
```
Usually we are interested in Lagrangian systems, where
```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) ,
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
abstract type IPRKMethod <: PRKMethod end

isexplicit(method::Union{IPRKMethod, Type{<:IPRKMethod}}) = false
isimplicit(method::Union{IPRKMethod, Type{<:IPRKMethod}}) = true

default_solver(::IPRKMethod) = Newton()
default_iguess(::IPRKMethod) = HermiteExtrapolation()


"""
Implicit Partitioned Runge-Kutta Method

```
IPRK(tableau)
```
"""
struct IPRK{TT <: PartitionedTableau} <: IPRKMethod
    tableau::TT
end

initmethod(method::IPRKMethod) = IPRK(method)

solversize(problem::AbstractProblemPODE, method::IPRK) = 2 * ndims(problem) * nstages(method)

function Base.show(io::IO, int::GeometricIntegrator{<:IPRK})
    print(io, "\nImplicit Partitioned Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end


@doc raw"""
Implicit partitioned Runge-Kutta integrator cache.

### Fields

* `x`: nonlinear solver solution vector
* `q̄`: solution at previous timestep
* `p̄`: momentum at previous timestep
* `Q`: internal stages of solution q
* `P`: internal stages of momentum p
* `V`: internal stages of vector field v = q̇
* `F`: internal stages of vector field f = ṗ
* `Y`: summed vector field of internal stages Q
* `Z`: summed vector field of internal stages P
"""
struct IPRKCache{ST,D,S,N} <: PODEIntegratorCache{ST,D}
    x::Vector{ST}

    q̄::Vector{ST}
    p̄::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}
    Z::Vector{Vector{ST}}

    function IPRKCache{ST,D,S,N}() where {ST,D,S,N}
        # create solver vector
        x = zeros(ST, N)

        # create previous solution vectors
        q̄ = zeros(ST, D)
        p̄ = zeros(ST, D)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        new(x, q̄, p̄, Q, P, V, F, Y, Z)
    end
end

function Cache{ST}(problem::EquationProblem, method::IPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IPRKCache{ST, D, S, solversize(problem, method)}(; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::IPRK) = IPRKCache{ST, ndims(problem), nstages(tableau(method)), solversize(problem, method)}

nlsolution(cache::IPRKCache) = cache.x

function reset!(cache::IPRKCache, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end


function initial_guess!(int::GeometricIntegrator{<:IPRK, <:AbstractProblemPODE})
    # get cache for nonlinear solution vector and internal stages
    local x = nlsolution(int)
    local Q = cache(int).Q
    local P = cache(int).P
    local V = cache(int).V
    local F = cache(int).F

    # compute initial guess for internal stages
    for i in eachstage(int)
        initialguess!(solstep(int).t̄ + timestep(int) * tableau(int).q.c[i], Q[i], P[i], V[i], F[i], solstep(int), problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(int)
        for k in 1:ndims(int)
            x[2*(ndims(int)*(i-1)+k-1)+1] = 0
            x[2*(ndims(int)*(i-1)+k-1)+2] = 0
            for j in eachstage(int)
                x[2*(ndims(int)*(i-1)+k-1)+1] += tableau(int).q.a[i,j] * V[j][k]
                x[2*(ndims(int)*(i-1)+k-1)+2] += tableau(int).p.a[i,j] * F[j][k]
            end
        end
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods from nonlinear solution vector
function components!(x::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemPODE}) where {ST}
    # get cache for internal stages
    local Q = cache(int, ST).Q
    local P = cache(int, ST).P
    local V = cache(int, ST).V
    local F = cache(int, ST).F
    local Y = cache(int, ST).Y
    local Z = cache(int, ST).Z
    local D = ndims(int)

    for i in eachstage(int)
        for k in 1:D
            # copy y to Y and Z
            Y[i][k] = x[2*(D*(i-1)+k-1)+1]
            Z[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            Q[i][k] = cache(int).q̄[k] + timestep(int) * Y[i][k]
            P[i][k] = cache(int).p̄[k] + timestep(int) * Z[i][k]
        end

        # compute v(Q,P) and f(Q,P)
        equations(int).v(V[i], solstep(int).t̄ + timestep(int) * tableau(int).q.c[i], Q[i], P[i], parameters(solstep(int)))
        equations(int).f(F[i], solstep(int).t̄ + timestep(int) * tableau(int).p.c[i], Q[i], P[i], parameters(solstep(int)))
    end
end


# Compute residual of implicit partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:IPRK, <:AbstractProblemPODE}) where {ST}
    # get cache for internal stages
    local V = cache(int, ST).V
    local F = cache(int, ST).F
    local Y = cache(int, ST).Y
    local Z = cache(int, ST).Z
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute b = - [(Y-AV), (Z-AF)]
    for i in eachindex(Y,Z)
        for k in eachindex(Y[i],Z[i])
            b[2*(D*(i-1)+k-1)+1] = - Y[i][k]
            b[2*(D*(i-1)+k-1)+2] = - Z[i][k]
            for j in eachindex(V,F)
                b[2*(D*(i-1)+k-1)+1] += tableau(int).q.a[i,j] * V[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(int).p.a[i,j] * F[j][k]
            end
        end
    end
end


function integrate_step!(int::GeometricIntegrator{<:IPRK, <:AbstractProblemPODE})
    # copy previous solution from solstep to cache
    reset!(cache(int), current(solstep(int))...)

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

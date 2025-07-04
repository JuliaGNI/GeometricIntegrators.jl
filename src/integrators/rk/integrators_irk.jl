@doc raw"""
Implicit Runge-Kutta method solving the explicit system of equations for [`ODE`](@ref)s
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
or the implicit systems of equations for [`IODE`](@ref)s and [`LODE`](@ref)s,
```math
\begin{aligned}
P_{n,i} &= \vartheta (Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , \\
F_{n,i} &= f (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} .
\end{aligned}
```
If `implicit_update` is set to `true`, the update is computed by solving
```math
\vartheta(q_{n+1}) = \vartheta(q_{n}) + h \sum \limits_{i=1}^{s} b_{i}  \, f (Q_{n,j}, V_{n,j}) ,
```
otherwise it is computed explicitly by
```math
q_{n+1} = q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
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
abstract type IRKMethod <: RKMethod end

isexplicit(method::Union{IRKMethod,Type{<:IRKMethod}}) = false
isimplicit(method::Union{IRKMethod,Type{<:IRKMethod}}) = true

implicit_update(::IRKMethod) = false

default_solver(::IRKMethod) = Newton()
default_iguess(::IRKMethod) = HermiteExtrapolation()

solversize(problem::AbstractProblemODE, method::IRKMethod) = ndims(problem) * nstages(method)


"""
Implicit Runge-Kutta Method

```
IRK(tableau)
```
"""
struct IRK{TT<:Tableau,ImplicitUpdate} <: IRKMethod
    tableau::TT

    function IRK(tableau::TT; implicit_update::Bool=false) where {TT<:Tableau}
        new{TT,implicit_update}(tableau)
    end
end

implicit_update(::IRK{TT,IU}) where {TT,IU} = IU

initmethod(method::IRKMethod) = IRK(method)


function Base.show(io::IO, int::GeometricIntegrator{<:IRK})
    print(io, "\nImplicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


@doc raw"""
Implicit Runge-Kutta integrator cache.

### Fields

* `x`: nonlinear solver solution vector
* `q̄`: solution at previous timestep
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Y`: vector field of internal stages
* `J`: Jacobi matrices for all internal stages
"""
struct IRKCache{DT,D,S} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    J::Vector{Matrix{DT}}

    function IRKCache{DT,D,S}() where {DT,D,S}
        x = zeros(DT, D * S)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        J = [zeros(DT, D, D) for _ in 1:S]
        new(x, Q, V, Y, J)
    end
end

nlsolution(cache::IRKCache) = cache.x

function Cache{ST}(problem::AbstractProblem, method::IRKMethod; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IRKCache{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblem, method::IRKMethod) = IRKCache{ST,ndims(problem),nstages(tableau(method))}


function internal_variables(method::IRKMethod, problem::AbstractProblemODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = ndims(problem)

    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    Y = create_internal_stage_vector(DT, D, S)

    # solver = get_solver_status(int.solver)

    (Q=Q, V=V, Y=Y)#, solver=solver)
end

function copy_internal_variables(solstep::SolutionStep, cache::IRKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :Y) && copyto!(internal(solstep).Y, cache.Y)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE})
    # compute initial guess for internal stages
    for i in eachstage(int)
        soltmp = (
            t=history.t[1] + timestep(int) * tableau(int).c[i],
            q=cache(int).Q[i],
            v=cache(int).V[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(int)
        offset = ndims(int) * (i - 1)
        for k in 1:ndims(int)
            cache(int).x[offset+k] = 0
            for j in eachstage(int)
                cache(int).x[offset+k] += timestep(int) * tableau(int).a[i, j] * cache(int).V[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)
    local D = ndims(int)

    # copy x to Y and compute Q = q + Y
    for i in eachindex(C.Q, C.Y)
        for k in eachindex(C.Q[i], C.Y[i])
            C.Y[i][k] = x[D*(i-1)+k]
            C.Q[i][k] = sol.q[k] + C.Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in eachindex(C.Q, C.V)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        equations(int).v(C.V[i], tᵢ, C.Q[i], params)
    end
end


# Compute stages of implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE}) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)
    local D = ndims(int)

    # compute b = - (Y-AV)
    for i in eachindex(C.Y)
        for k in eachindex(C.Y[i])
            y1 = y2 = zero(ST)
            for j in eachindex(C.V)
                y1 += tableau(int).a[i, j] * C.V[j][k]
                y2 += tableau(int).â[i, j] * C.V[j][k]
            end
            b[D*(i-1)+k] = C.Y[i][k] - timestep(int) * (y1 + y2)
        end
    end
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages of implicit Runge-Kutta methods from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute right-hand side b of nonlinear solver
    residual!(b, sol, params, int)
end


function update!(sol, params, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE}, DT)
    # compute final update
    update!(sol.q, cache(int, DT).V, tableau(int), timestep(int))
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:IRK,<:AbstractProblemODE})
    # call nonlinear solver
    solve!(solver(int), nlsolution(int), (sol, params, int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end

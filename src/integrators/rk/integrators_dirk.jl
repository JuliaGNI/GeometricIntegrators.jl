@doc raw"""
Diagonally implicit Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
"""
abstract type DIRKMethod <: IRKMethod end

default_solver(::DIRKMethod) = Newton()
default_iguess(::DIRKMethod) = HermiteExtrapolation()


"""
Diagonally Implicit Runge-Kutta Method

```
DIRK(tableau)
```
"""
struct DIRK{TT<:Tableau} <: DIRKMethod
    tableau::TT
end

initmethod(method::DIRKMethod) = DIRK(method)

solversize(problem::AbstractProblemODE, ::DIRKMethod) = ndims(problem)


function Base.show(io::IO, int::GeometricIntegrator{<:DIRK})
    print(io, "\nDiagonally Implicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


struct SingleStageSolvers{ST} <: NonlinearSolver
    solvers::ST
    SingleStageSolvers(solvers...) = new{typeof(solvers)}(solvers)
end

Base.getindex(s::SingleStageSolvers, args...) = getindex(s.solvers, args...)


function initsolver(::NewtonMethod, method::DIRK, caches::CacheDict; kwargs...)
    SingleStageSolvers([NewtonSolver(zero(cache(caches).x[i]), residual!, zero(cache(caches).x[i]); linesearch=Backtracking(), kwargs...) for i in eachstage(method)]...)
end


"""
Diagonally implicit Runge-Kutta integrator cache.

### Fields

* `x`: nonlinear solver solution vector
* `Q`: internal stages of solution q
* `V`: internal stages of vector field v = q̇
* `Y`: summed vector field of internal stages Q
"""
struct DIRKCache{DT,D,S} <: ODEIntegratorCache{DT,D}
    x::Vector{Vector{DT}}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function DIRKCache{DT,D,S}() where {DT,D,S}
        x = create_internal_stage_vector(DT, D, S)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        new(x, Q, V, Y)
    end
end

function Cache{ST}(problem::EquationProblem, method::DIRKMethod; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    DIRKCache{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::DIRKMethod) = DIRKCache{ST,ndims(problem),nstages(tableau(method))}

nlsolution(cache::DIRKCache, i) = cache.x[i]


function internal_variables(method::DIRKMethod, problem::AbstractProblemODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = ndims(problem)

    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    Y = create_internal_stage_vector(DT, D, S)

    # solver = get_solver_status(int.solver)

    (Q=Q, V=V, Y=Y)#, solver=solver)
end

function copy_internal_variables!(solstep::SolutionStep, cache::DIRKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :Y) && copyto!(internal(solstep).Y, cache.Y)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:DIRK})
    for i in eachstage(int)
        soltmp = (
            t=history.t[1] + timestep(int) * tableau(int).c[i],
            q=cache(int).Q[i],
            q̇=cache(int).V[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end

    for i in eachstage(int)
        for k in 1:ndims(int)
            cache(int).x[i][k] = 0
            for j in eachstage(int)
                cache(int).x[i][k] += timestep(int) * tableau(int).a[i, j] * cache(int).V[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DIRK}, i) where {ST}
    # get cache for internal stages
    local C = cache(int, ST)

    # copy x to Y
    for k in eachindex(C.Y[i])
        C.Y[i][k] = x[k]
    end

    # compute Q = q + Y
    C.Q[i] .= sol.q .+ C.Y[i]

    # compute V = v(Q)
    equations(int).v(C.V[i], sol.t + timestep(int) * (tableau(int).c[i] - 1), C.Q[i], params)
end


function residual!(b::AbstractVector{ST}, int::GeometricIntegrator{<:DIRK}, i) where {ST}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    local C = cache(int, ST)

    # compute b = - (Y-AV)
    for k in 1:ndims(int)
        y1 = tableau(int).a[i, i] * C.V[i][k]
        y2 = tableau(int).â[i, i] * C.V[i][k]
        for j in 1:i-1
            y1 += tableau(int).a[i, j] * C.V[j][k]
            y2 += tableau(int).â[i, j] * C.V[j][k]
        end
        b[k] = C.Y[i][k] - timestep(int) * (y1 + y2)
    end
end


# Compute stages of diagonally implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DIRK}, i) where {ST}
    # compute stages from nonlinear solver solution x
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int, i)

    # compute residual vector
    residual!(b, int, i)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:DIRK,<:AbstractProblemODE})
    # copy previous solution from solstep to cache
    reset!(cache(int), sol...)

    # consecutively solve for all stages
    for i in eachstage(int)
        # call nonlinear solver
        solve!(solver(int)[i], nlsolution(cache(int), i), (sol, params, int, i))

        # print solver status
        # println(status(solvers[i]))

        # check if solution contains NaNs or error bounds are violated
        # println(meets_stopping_criteria(status(solvers[i])))

        # compute vector field at internal stages
        components!(nlsolution(cache(int), i), sol, params, int, i)
    end

    # compute final update
    update!(sol.q, cache(int).V, tableau(int), timestep(int))
end

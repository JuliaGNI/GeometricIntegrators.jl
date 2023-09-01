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
struct DIRK{TT <: Tableau} <: DIRKMethod
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
    SingleStageSolvers([NewtonSolver(zero(cache(caches).x[i]), zero(cache(caches).x[i]); linesearch = Backtracking(), config = Options(min_iterations = 1)) for i in eachstage(method)]...)
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

@inline CacheType(ST, problem::EquationProblem, method::DIRKMethod) = DIRKCache{ST, ndims(problem), nstages(tableau(method))}

nlsolution(cache::DIRKCache, i) = cache.x[i]


function initial_guess!(int::GeometricIntegrator{<:DIRK})
    for i in eachstage(int)
        initialguess!(solstep(int).t̄ + timestep(int) * tableau(int).c[i], cache(int).Q[i], cache(int).V[i], solstep(int), problem(int), iguess(int))
    end

    for i in eachstage(int)
        for k in 1:ndims(int)
            cache(int).x[i][k] = 0
            for j in eachstage(int)
                cache(int).x[i][k] += timestep(int) * tableau(int).a[i,j] * cache(int).V[j][k]
            end
        end
    end
end


function components!(x::AbstractVector{ST}, int::GeometricIntegrator{<:DIRK}, i) where {ST}
    # time of i-th stage
    local tᵢ::timetype(problem(int))

    # get cache for internal stages
    local Q::Vector{Vector{ST}} = cache(int, ST).Q
    local V::Vector{Vector{ST}} = cache(int, ST).V
    local Y::Vector{Vector{ST}} = cache(int, ST).Y

    # copy x to Y and compute Q = q + Y
    for k in 1:ndims(int)
        Y[i][k] = x[k]
        Q[i][k] = solstep(int).q[k] + Y[i][k]
    end

    # compute V = v(Q)
    tᵢ = solstep(int).t̄ + timestep(int) * tableau(int).c[i]
    equations(int).v(V[i], tᵢ, Q[i])
end


function residual!(b::AbstractVector{ST}, int::GeometricIntegrator{<:DIRK}, i) where {ST}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    local V::Vector{Vector{ST}} = cache(int, ST).V
    local Y::Vector{Vector{ST}} = cache(int, ST).Y

    # compute b = - (Y-AV)
    for k in 1:ndims(int)
        y1 = tableau(int).a[i,i] * V[i][k]
        y2 = tableau(int).â[i,i] * V[i][k]
        for j in 1:i-1
            y1 += tableau(int).a[i,j] * V[j][k]
            y2 += tableau(int).â[i,j] * V[j][k]
        end
        b[k] = - Y[i][k] + timestep(int) * (y1 + y2)
    end
end


# Compute stages of diagonally implicit Runge-Kutta methods.
function residual!(b, x, int::GeometricIntegrator{<:DIRK}, i)
    # compute stages from nonlinear solver solution x
    components!(x, int, i)

    # compute residual vector
    residual!(b, int, i)
end


function integrate_step!(int::GeometricIntegrator{<:DIRK, <:AbstractProblemODE})
    # consecutively solve for all stages
    for i in eachstage(int)
        # call nonlinear solver
        solve!(nlsolution(cache(int), i), (b,x) -> residual!(b, x, int, i), solver(int)[i])

        # print solver status
        # println(status(solvers[i]))

        # check if solution contains NaNs or error bounds are violated
        # println(meets_stopping_criteria(status(solvers[i])))

        # compute vector field at internal stages
        components!(nlsolution(cache(int), i), int, i)
    end

    # compute final update
    update!(solstep(int), cache(int).V, tableau(int), timestep(int))
end

@doc raw"""
Implicit Runge-Kutta integrator cache.

### Fields

* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Y`: vector field of internal stages
* `J`: Jacobi matrices for all internal stages
"""
struct IntegratorCacheIRK{DT,D,S} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    J::Vector{Matrix{DT}}

    function IntegratorCacheIRK{DT,D,S}() where {DT,D,S}
        x = zeros(DT, D*S)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        J = [zeros(DT,D,D) for _ in 1:S]
        new(x, Q, V, Y, J)
    end
end

function Cache{ST}(problem::GeometricProblem, method::IRKMethod; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheIRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::GeometricProblem, method::IRKMethod) = IntegratorCacheIRK{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Implicit Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} .
\end{aligned}
```
"""
const IntegratorIRK{DT,TT} = Integrator{<:ODEProblem{DT,TT}, <:IRKMethod}


initmethod(method::IRK) = method
initmethod(method::IRKMethod) = IRK(method)

default_solver(::IRKMethod) = Newton()
default_iguess(::IRKMethod) = HermiteExtrapolation()


function initsolver(::Newton, solstep::SolutionStepODE{DT}, problem::ODEProblem, method::IRKMethod, caches::CacheDict) where {DT}
    # create wrapper function f!(b,x)
    f! = (b,x) -> function_stages!(b, x, solstep, problem, method, caches)

    # create nonlinear solver
    NewtonSolver(zero(caches[DT].x), zero(caches[DT].x), f!; linesearch = Backtracking(), config = Options(min_iterations = 1))
end


function Base.show(io::IO, int::IntegratorIRK)
    print(io, "\nImplicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


function initial_guess!(
    solstep::SolutionStepODE{DT}, 
    problem::ODEProblem,
    method::IRKMethod, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}
    
    # obtain cache
    cache = caches[DT]

    local offset::Int

    # compute initial guess for internal stages
    for i in eachstage(method)
        initialguess!(solstep.t̄[1] + timestep(problem) * tableau(method).c[i], cache.Q[i], cache.V[i], solstep, problem, iguess)
    end

    # assemble initial guess for nonlinear solver solution vector
    for i in eachstage(method)
        offset = ndims(problem)*(i-1)
        for k in 1:ndims(problem)
            cache.x[offset+k] = 0
            for j in eachstage(method)
                cache.x[offset+k] += timestep(problem) * tableau(method).a[i,j] * cache.V[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, solstep::SolutionStepODE{DT,TT}, problem::ODEProblem, method::IRKMethod, caches::CacheDict) where {ST,DT,TT}
    local tᵢ::TT

    local Q = caches[ST].Q
    local V = caches[ST].V
    local Y = caches[ST].Y
    local D = ndims(problem)

    # copy x to Y and compute Q = q + Y
    for i in eachindex(Q,Y)
        for k in eachindex(Q[i],Y[i])
            Y[i][k] = x[D*(i-1)+k]
            Q[i][k] = solstep.q̄[1][k] + Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in eachindex(Q,V)
        tᵢ = solstep.t̄[1] + timestep(problem) * tableau(method).c[i]
        functions(problem).v(V[i], tᵢ, Q[i])
    end
end

# Compute stages of implicit Runge-Kutta methods.
function function_stages!(b::Vector{ST}, x::Vector{ST}, solstep::SolutionStepODE, problem::ODEProblem, method::IRKMethod, caches::CacheDict) where {ST}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    local V = caches[ST].V
    local Y = caches[ST].Y
    local D = ndims(problem)

    # compute stages from nonlinear solver solution x
    compute_stages!(x, solstep, problem, method, caches)

    # compute b = - (Y-AV)
    for i in eachindex(Y)
        for k in eachindex(Y[i])
            y1 = y2 = 0
            for j in eachindex(V)
                y1 += tableau(method).a[i,j] * V[j][k]
                y2 += tableau(method).â[i,j] * V[j][k]
            end
            b[D*(i-1)+k] = - Y[i][k] + timestep(problem) * (y1 + y2)
        end
    end
end


function integrate_step!(
    solstep::SolutionStepODE{DT,TT},
    problem::ODEProblem{DT,TT}, 
    method::IRKMethod,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute vector field at internal stages
    compute_stages!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    update_solution!(solstep.q, solstep.q̄[1], solstep.q̃, caches[DT].V, tableau(method).b, tableau(method).b̂, timestep(problem))

    # update vector field for initial guess
    update_vector_fields!(solstep, problem)
end

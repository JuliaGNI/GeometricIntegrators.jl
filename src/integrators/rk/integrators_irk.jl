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
    q̄::Vector{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    J::Vector{Matrix{DT}}

    function IntegratorCacheIRK{DT,D,S}() where {DT,D,S}
        x = zeros(DT, D*S)
        q̄ = zeros(DT, D)
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        Y = create_internal_stage_vector(DT, D, S)
        J = [zeros(DT,D,D) for _ in 1:S]
        new(x, q̄, Q, V, Y, J)
    end
end

nlsolution(cache::IntegratorCacheIRK) = cache.x
reset!(cache::IntegratorCacheIRK, t, q, λ = missing) = copyto!(cache.q̄, q)


function Cache{ST}(problem::AbstractProblem, method::IRKMethod; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheIRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblem, method::IRKMethod) = IntegratorCacheIRK{ST, ndims(problem), nstages(tableau(method))}


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
const IntegratorIRK{DT,TT} = Integrator{<:Union{ODEProblem{DT,TT}, DAEProblem{DT,TT}, SubstepProblem{DT,TT}}, <:IRKMethod}

solversize(problem::Union{ODEProblem, DAEProblem, SubstepProblem}, method::IRKMethod) =
    ndims(problem) * nstages(method)

initmethod(method::IRK) = method
initmethod(method::IRKMethod) = IRK(method)

default_solver(::IRKMethod) = Newton()
default_iguess(::IRKMethod) = HermiteExtrapolation()


function Base.show(io::IO, int::IntegratorIRK)
    print(io, "\nImplicit Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
    # print(io, reference(int.params.tab))
end


function initial_guess!(
    solstep::Union{SolutionStepODE{DT}, SolutionStepDAE{DT}, SubstepProblem{DT}},
    problem::Union{ODEProblem, DAEProblem, SubstepProblem},
    method::IRKMethod, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}
    
    # obtain cache
    cache = caches[DT]

    local offset::Int

    # compute initial guess for internal stages
    for i in eachstage(method)
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).c[i], cache.Q[i], cache.V[i], solstep, problem, iguess)
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


function components!(x::AbstractVector{ST}, int::IntegratorIRK) where {ST}
    local tᵢ::timetype(problem(int))

    local q̄ = cache(int, ST).q̄
    local Q = cache(int, ST).Q
    local V = cache(int, ST).V
    local Y = cache(int, ST).Y
    local D = ndims(int)

    # copy x to Y and compute Q = q + Y
    for i in eachindex(Q,Y)
        for k in eachindex(Q[i],Y[i])
            Y[i][k] = x[D*(i-1)+k]
            Q[i][k] = q̄[k] + Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in eachindex(Q,V)
        tᵢ = solstep(int).t̄ + timestep(int) * tableau(int).c[i]
        equations(int).v(V[i], tᵢ, Q[i])
    end
end


function residual!(b::AbstractVector{ST}, int::IntegratorIRK) where {ST}
    # temporary variables
    local y1::ST
    local y2::ST

    # get cache for internal stages
    local V = cache(int, ST).V
    local Y = cache(int, ST).Y
    local D = ndims(int)

    # compute b = - (Y-AV)
    for i in eachindex(Y)
        for k in eachindex(Y[i])
            y1 = y2 = 0
            for j in eachindex(V)
                y1 += tableau(int).a[i,j] * V[j][k]
                y2 += tableau(int).â[i,j] * V[j][k]
            end
            b[D*(i-1)+k] = - Y[i][k] + timestep(int) * (y1 + y2)
        end
    end
end


# Compute stages of implicit Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::IntegratorIRK) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::IntegratorIRK) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    update!(solstep(int), cache(int, DT).V, tableau(int), timestep(int))
end


function integrate_step!(int::IntegratorIRK)
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(nlsolution(int), int)
end

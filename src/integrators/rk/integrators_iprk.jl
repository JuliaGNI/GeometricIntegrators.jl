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
struct IntegratorCacheIPRK{ST,D,S} <: PODEIntegratorCache{ST,D}
    x::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}
    Z::Vector{Vector{ST}}

    function IntegratorCacheIPRK{ST,D,S}() where {ST,D,S}
        # create solver vector
        x = zeros(ST, 2*D*S)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        new(x, Q, P, V, F, Y, Z)
    end
end

nlsolution(cache::IntegratorCacheIPRK) = cache.x


function Cache{ST}(problem::GeometricProblem, method::IPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheIPRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::GeometricProblem, method::IPRK) = IntegratorCacheIPRK{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Implicit partitioned Runge-Kutta integrator solving the system
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
"""
const IntegratorIPRK{DT,TT} = Integrator{<:Union{PODEProblem{DT,TT},HODEProblem{DT,TT}}, <:IPRK}

solversize(problem::Union{PODEProblem,HODEProblem}, method::IPRK) =
    2 * ndims(problem) * nstages(method)

initmethod(method::IPRK) = method
initmethod(method::IPRKMethod) = IPRK(method)

default_solver(::IPRKMethod) = Newton()
default_iguess(::IPRKMethod) = HermiteExtrapolation()


function Base.show(io::IO, int::IntegratorIPRK)
    print(io, "\nImplicit Partitioned Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end


function initial_guess!(
    solstep::SolutionStepPODE{DT}, 
    problem::Union{PODEProblem,HODEProblem},
    method::IPRKMethod, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}

    local Q = caches[DT].Q
    local P = caches[DT].P
    local V = caches[DT].V
    local F = caches[DT].F

    for i in eachstage(method)
        initialguess!(solstep.t̄ + timestep(problem) * tableau(method).q.c[i], Q[i], P[i], V[i], F[i], solstep, problem, iguess)
    end
    for i in eachstage(method)
        for k in 1:ndims(problem)
            caches[DT].x[2*(ndims(problem)*(i-1)+k-1)+1] = 0
            caches[DT].x[2*(ndims(problem)*(i-1)+k-1)+2] = 0
            for j in eachstage(method)
                caches[DT].x[2*(ndims(problem)*(i-1)+k-1)+1] += tableau(method).q.a[i,j] * V[j][k]
                caches[DT].x[2*(ndims(problem)*(i-1)+k-1)+2] += tableau(method).p.a[i,j] * F[j][k]
            end
        end
    end
end


function components!(
    x::AbstractVector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{PODEProblem,HODEProblem},
    method::IPRKMethod,
    caches::CacheDict) where {ST,DT,TT}

    local Q = caches[ST].Q
    local P = caches[ST].P
    local V = caches[ST].V
    local F = caches[ST].F
    local Y = caches[ST].Y
    local Z = caches[ST].Z

    local tqᵢ::TT
    local tpᵢ::TT

    local D = ndims(problem)

    for i in eachstage(method)
        for k in 1:D
            # copy y to Y and Z
            Y[i][k] = x[2*(D*(i-1)+k-1)+1]
            Z[i][k] = x[2*(D*(i-1)+k-1)+2]

            # compute Q and P
            Q[i][k] = solstep.q̄[1][k] + timestep(problem) * Y[i][k]
            P[i][k] = solstep.p̄[1][k] + timestep(problem) * Z[i][k]
        end

        # compute time of internal stage
        tqᵢ = solstep.t̄ + timestep(problem) * tableau(method).q.c[i]
        tpᵢ = solstep.t̄ + timestep(problem) * tableau(method).p.c[i]

        # compute v(Q,P) and f(Q,P)
        functions(problem).v(V[i], tqᵢ, Q[i], P[i])
        functions(problem).f(F[i], tpᵢ, Q[i], P[i])
    end
end


# Compute stages of implicit partitioned Runge-Kutta methods.
function residual!(
    b::AbstractVector{ST},
    x::AbstractVector{ST},
    solstep::SolutionStepPODE,
    problem::Union{PODEProblem,HODEProblem},
    method::IPRKMethod,
    caches::CacheDict) where {ST}

    # get cache for internal stages
    local V = caches[ST].V
    local F = caches[ST].F
    local Y = caches[ST].Y
    local Z = caches[ST].Z

    local D = ndims(problem)

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)

    # compute b = - [(Y-AV), (Z-AF)]
    for i in eachstage(method)
        for k in 1:D
            b[2*(D*(i-1)+k-1)+1] = - Y[i][k]
            b[2*(D*(i-1)+k-1)+2] = - Z[i][k]
            for j in eachstage(method)
                b[2*(D*(i-1)+k-1)+1] += tableau(method).q.a[i,j] * V[j][k]
                b[2*(D*(i-1)+k-1)+2] += tableau(method).p.a[i,j] * F[j][k]
            end
        end
    end
end


function integrate_step!(
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{PODEProblem{DT,TT},HODEProblem{DT,TT}},
    method::IPRK,
    caches::CacheDict,
    solver::NonlinearSolver) where {DT,TT}

    # call nonlinear solver
    solve!(caches[DT].x, (b,x) -> residual!(b, x, solstep, problem, method, caches), solver)

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute vector fields at internal stages
    components!(caches[DT].x, solstep, problem, method, caches)

    # compute final update
    update!(solstep, caches[DT].V, caches[DT].F, tableau(method), timestep(problem))
end

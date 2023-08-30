@doc raw"""
Trapezoidal Variational Integrator in position-momentum form.

We consider a discrete Lagrangian of the form
```math
L_d (q_{n}, q_{n+1}) = \frac{h}{2} \bigg[ L \bigg( q_{n}, \frac{q_{n+1} - q_{n}}{h} \bigg) + L \bigg( q_{n+1}, \frac{q_{n+1} - q_{n}}{h} \bigg) \bigg] ,
```
where $q_{n}$ approximates the solution $q(t_{n})$.
The Euler-Lagrange equations are computed as:
```math
\begin{aligned}
p_{n  } &=           -  D_{1} L_{d} (q_{n}, q_{n+1}) = \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( q_{n}, \frac{q_{n+1} - q_{n}}{h} \bigg) + \frac{1}{2} \frac{\partial L}{\partial v} \bigg( q_{n}, \frac{q_{n+1} - q_{n}}{h} \bigg) + \frac{1}{2} \frac{\partial L}{\partial v} \bigg( q_{n+1}, \frac{q_{n+1} - q_{n}}{h} \bigg) , \\
p_{n+1} &= \hphantom{-} D_{2} L_{d} (q_{n}, q_{n+1}) = \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( q_{n}, \frac{q_{n+1} - q_{n}}{h} \bigg) + \frac{1}{2} \frac{\partial L}{\partial v} \bigg( q_{n}, \frac{q_{n+1} - q_{n}}{h} \bigg) + \frac{1}{2} \frac{\partial L}{\partial v} \bigg( q_{n+1}, \frac{q_{n+1} - q_{n}}{h} \bigg) .
\end{aligned}
```
The first equation can be solved implicitly for $q_{n+1}$ given $(q_{n}, p_{n})$.
The second equation can be used to explicitly compute $p_{n+1}$.
"""
struct PMVItrapezoidal <: VIMethod end

isexplicit(method::PMVItrapezoidal) = false
isimplicit(method::PMVItrapezoidal) = true
issymmetric(method::PMVItrapezoidal) = true
issymplectic(method::PMVItrapezoidal) = true


const PMVItrapezoidalIntegrator{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:PMVItrapezoidal}

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::PMVItrapezoidal; kwargs...) where {ST}
    IntegratorCachePMVI{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, ::PMVItrapezoidal) = IntegratorCachePMVI{ST, ndims(problem)}

function Base.show(io::IO, int::PMVItrapezoidalIntegrator)
    print(io, "\nTrapezoidal variational integrator in position-momentum form with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::PMVItrapezoidal,
    caches::CacheDict) where {ST,DT,TT}

    # get cache and dimension
    cache = caches[ST]
    D = ndims(problem)

    # set some local variables for convenience and clarity
    local t̄ = solstep.t̄
    local t = solstep.t̄ + timestep(problem)
    
    # copy x to q
    cache.q .= x[1:D]

    # compute v
    cache.ṽ .= (cache.q .- cache.q̄) ./ timestep(problem)
 
    # compute Θ = ϑ(q,ṽ) and f = f(q,ṽ)
    functions(problem).ϑ(cache.θ̄, t̄, cache.q̄, cache.ṽ)
    functions(problem).ϑ(cache.θ, t, cache.q, cache.ṽ)
    functions(problem).f(cache.f̄, t̄, cache.q̄, cache.ṽ)
    functions(problem).f(cache.f, t, cache.q, cache.ṽ)

    # compute p
    cache.θ̃ .= (cache.θ .+ cache.θ̄) ./ 2
    cache.p .= cache.p̄ .+ timestep(problem) .* (cache.f .+ cache.f̄) ./ 2
end


function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::PMVItrapezoidal,
    caches::CacheDict) where {ST}

    # get cache for internal stages
    cache = caches[ST]

    # copy previous solution from solstep to cache
    reset!(cache, current(solstep)...)

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)

    # compute b
    b .= cache.θ̃ .- cache.p̄ .- timestep(problem) .* cache.f̄ ./ 2
end

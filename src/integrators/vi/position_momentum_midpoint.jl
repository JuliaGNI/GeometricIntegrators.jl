@doc raw"""
Midpoint Variational Integrator in position-momentum form.

We consider a discrete Lagrangian of the form
```math
L_d (q_{n}, q_{n+1}) = h \, L \bigg( \frac{q_{n} + q_{n+1}}{2}, \frac{q_{n+1} - q_{n}}{h} \bigg) ,
```
where $q_{n}$ approximates the solution $q(t_{n})$.
The Euler-Lagrange equations are computed as:
```math
\begin{aligned}
p_{n  } &=           -  D_{1} L_{d} (q_{n}, q_{n+1}) = \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, \frac{q_{n+1} - q_{n}}{h} \bigg) + \frac{\partial L}{\partial v} \bigg( \frac{q_{n} + q_{n+1}}{2}, \frac{q_{n+1} - q_{n}}{h} \bigg) , \\
p_{n+1} &= \hphantom{-} D_{2} L_{d} (q_{n}, q_{n+1}) = \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, \frac{q_{n+1} - q_{n}}{h} \bigg) + \frac{\partial L}{\partial v} \bigg( \frac{q_{n} + q_{n+1}}{2}, \frac{q_{n+1} - q_{n}}{h} \bigg) .
\end{aligned}
```
The first equation can be solved implicitly for $q_{n+1}$ given $(q_{n}, p_{n})$.
The second equation can be used to explicitly compute $p_{n+1}$.    
"""
struct PMVImidpoint <: VIMethod end

isexplicit(method::PMVImidpoint) = false
isimplicit(method::PMVImidpoint) = true
issymmetric(method::PMVImidpoint) = true
issymplectic(method::PMVImidpoint) = true


const PMVImidpointIntegrator{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:PMVImidpoint}

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::PMVImidpoint; kwargs...) where {ST}
    IntegratorCachePMVI{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, ::PMVImidpoint) = IntegratorCachePMVI{ST, ndims(problem)}

function Base.show(io::IO, int::PMVImidpointIntegrator)
    print(io, "\nMidpoint variational integrator in position-momentum form with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::PMVImidpoint,
    caches::CacheDict) where {ST,DT,TT}

    # get cache and dimension
    cache = caches[ST]
    D = ndims(problem)

    # set some local variables for convenience and clarity
    local t̃ = solstep.t̄ + timestep(problem) / 2
    
    # copy x to q
    cache.q .= x[1:D]

    # compute q̃ and ṽ
    cache.q̃ .= (cache.q .+ cache.q̄) ./ 2
    cache.ṽ .= (cache.q .- cache.q̄) ./ timestep(problem)
 
    # compute Θ̃ = ϑ(q̃,ṽ) and f̃ = f(q̃,ṽ)
    functions(problem).ϑ(cache.θ̃, t̃, cache.q̃, cache.ṽ)
    functions(problem).f(cache.f̃, t̃, cache.q̃, cache.ṽ)

    # compute p
    cache.p .= cache.θ̃ .+ timestep(problem) .* cache.f̃ ./ 2
end


function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::PMVImidpoint,
    caches::CacheDict) where {ST}

    # get cache for internal stages
    cache = caches[ST]

    # copy previous solution from solstep to cache
    reset!(cache, current(solstep)...)

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)

    # compute b
    b .= cache.θ̃ .- cache.p̄ .- timestep(problem) .* cache.f̃ ./ 2
end

@doc raw"""
Hamilton-Pontryagin Integrator using midpoint quadrature.

We consider a discrete Lagrangian of the form
```math
L_d (q_{n}, q_{n+1}) = h \, L \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) ,
```
where $q_{n}$ approximates the solution $q(t_n)$ and $v_{n+1/2}$ is the velocity, which is assumed to be constant in the interval $[t_{n}, t_{n+1}]$.
The discrete Hamilton-Pontryagin action reads
```math
A_d [q_d] = \sum \limits_{n=0}^{N-1} \bigg[ L_{d}^{\alpha,\beta} (q_{n}, q_{n+1}, v_{n+1/2}) + \left< p_{n+1/2} , \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) - v_{n+1/2} \right> \bigg] ,
```
where $\phi_h$ is a map that computes the velocity $v_{n+1/2}$ as a function of $q_{n}$, $q_{n+1}$ and a set of parameters $a_{n,n+1}$.
A trivial example of such a map that does not depend on any parameters $a_{n,n+1}$ is
```math
\phi_h (q_{n}, q_{n+1}; a_{n,n+1}) = \frac{q_{n+1} - q_{n}}{h} .
```
In order to solve the discrete Euler-Lagrange equations, the user needs to specify the map $\phi_h$, its gradients with respect to $q_{n}$ and $q_{n+1}$, denoted by $D_1 \phi_h$ and $D_2 \phi_h$, respectively, the gradient with respect to the parameters, denoted by $D_a \phi_h$, and an initial set of parameters $a_{0}$.

The equations of motion, that are solved by this integrator, is computed as:
```math
\begin{aligned}
0 &= \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n-1} + q_{n}}{2}, v_{n-1/2} \bigg)
   + \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) \\
  &+ h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2}
   + h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n-1/2} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
p_{n+1/2} &= \frac{\partial L}{\partial v} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) , \\
v_{n+1/2} &= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) .
\end{aligned}
```
Upon defining the momentum
```math
p_{n}
= h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n-1/2}
+ \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n-1} + q_{n}}{2}, v_{n-1/2} \bigg) ,
```
we can rewrite the equations of motion as
```math
\begin{aligned}
0 &= p_{n}
   + \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg)
   + h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
v_{n+1/2} &= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) , \\
p_{n+1/2} &= \frac{\partial L}{\partial v} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) , \\
p_{n+1}
&= h \, D_2 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2}
 + \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) .
\end{aligned}
```
Given $(q_{n}, p_{n})$, the first four equations are solved for $q_{n+1}$, where $v_{n+1/2}$, $p_{n+1/2}$, and $a_{n,n+1}$ are treated as internal variables similar to the internal stages of a Runge-Kutta method.
The last equation provides an explicit update for $p_{n+1}$.
"""
struct HPImidpoint{ϕT, D₁ϕT, D₂ϕT, DₐϕT, PT} <: HPIMethod
    ϕ::ϕT
    D₁ϕ::D₁ϕT
    D₂ϕ::D₂ϕT
    Dₐϕ::DₐϕT
    params::PT
end

nparams(method::HPImidpoint) = length(method.params)

isexplicit(method::HPImidpoint) = false
isimplicit(method::HPImidpoint) = true
issymmetric(method::HPImidpoint) = missing
issymplectic(method::HPImidpoint) = true


const HPImidpointIntegrator{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:HPImidpoint}

function Cache{ST}(problem::Union{IODEProblem,LODEProblem}, method::HPImidpoint; kwargs...) where {ST}
    IntegratorCacheHPI{ST, ndims(problem), nparams(method)}(; kwargs...)
end

@inline CacheType(ST, problem::Union{IODEProblem,LODEProblem}, method::HPImidpoint) = IntegratorCacheHPI{ST, ndims(problem), nparams(method)}

function Base.show(io::IO, int::HPImidpointIntegrator)
    print(io, "\nHamilton-Pontryagin Integrator using midpoint quadrature with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function components!(
    x::Vector{ST},
    solstep::SolutionStepPODE{DT,TT},
    problem::Union{IODEProblem,LODEProblem},
    method::HPImidpoint,
    caches::CacheDict) where {ST,DT,TT}

    # get cache and dimension
    cache = caches[ST]
    D = ndims(problem)
    A = nparams(method)

    # set some local variables for convenience and clarity
    local t̃ = solstep.t̄ + timestep(problem) / 2
    
    # copy x to q
    cache.q .= x[1:D]
    cache.a .= x[D+1:D+A]

    # compute q̃
    cache.q̃ .= (cache.q .+ cache.q̄) ./ 2

    # compute v
    method.ϕ(cache.ṽ, cache.q̄, cache.q, cache.a, timestep(problem))
 
    # compute Θ̃ = ϑ(q̃,ṽ) and f̃ = f(q̃,ṽ)
    functions(problem).ϑ(cache.θ̃, t̃, cache.q̃, cache.ṽ)
    functions(problem).f(cache.f̃, t̃, cache.q̃, cache.ṽ)

    # compute derivatives of ϕ
    method.D₁ϕ(cache.D₁ϕ, cache.q̄, cache.q, cache.a, timestep(problem))
    method.D₂ϕ(cache.D₂ϕ, cache.q̄, cache.q, cache.a, timestep(problem))
    method.Dₐϕ(cache.Dₐϕ, cache.q̄, cache.q, cache.a, timestep(problem))

    # compute p
    cache.p .= timestep(problem) .* cache.f̃ ./ 2
    for i in 1:D
        for j in 1:D
            cache.p[i] += timestep(problem) * cache.D₂ϕ[i,j] * cache.θ̃[j]
        end
    end
end


function residual!(
    b::Vector{ST},
    x::Vector{ST},
    solstep::SolutionStepPODE,
    problem::Union{IODEProblem,LODEProblem},
    method::HPImidpoint,
    caches::CacheDict) where {ST}

    # get cache for internal stages
    cache = caches[ST]
    D = ndims(problem)
    A = nparams(method)

    # copy previous solution from solstep to cache
    reset!(cache, current(solstep)...)

    # compute stages from nonlinear solver solution x
    components!(x, solstep, problem, method, caches)

    # compute b
    for i in 1:D
        b[i] = cache.p̄[i] + timestep(problem) * cache.f̃[i] / 2
        for j in 1:D
            b[i] += timestep(problem) * cache.D₁ϕ[i,j] * cache.θ̃[j]
        end
    end
    for i in 1:A
        b[D+i] = 0
        for j in 1:D
            b[D+i] += cache.Dₐϕ[i,j] * cache.θ̃[j]
        end
    end
end
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
struct PMVImidpoint <: PMVIMethod end

issymmetric(method::PMVImidpoint) = true


function Base.show(io::IO, int::GeometricIntegrator{<:PMVImidpoint})
    print(io, "\nMidpoint variational integrator in position-momentum form with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function components!(x::Vector{ST}, int::GeometricIntegrator{<:PMVImidpoint}) where {ST}
    # set some local variables for convenience and clarity
    local t̃ = solstep(int).t̄ + timestep(int) / 2
    
    # copy x to q
    cache(int, ST).q .= x[1:ndims(int)]

    # compute q̃ and ṽ
    cache(int, ST).q̃ .= (cache(int, ST).q .+ cache(int, ST).q̄) ./ 2
    cache(int, ST).ṽ .= (cache(int, ST).q .- cache(int, ST).q̄) ./ timestep(int)
 
    # compute Θ̃ = ϑ(q̃,ṽ) and f̃ = f(q̃,ṽ)
    equations(int).ϑ(cache(int, ST).θ̃, t̃, cache(int, ST).q̃, cache(int, ST).ṽ)
    equations(int).f(cache(int, ST).f̃, t̃, cache(int, ST).q̃, cache(int, ST).ṽ)

    # compute p
    cache(int, ST).p .= cache(int, ST).p̄ .+ timestep(int) .* cache(int, ST).f̃
end


function residual!(b::Vector{ST}, int::GeometricIntegrator{<:PMVImidpoint}) where {ST}
    # compute b
    b .= cache(int, ST).θ̃ .- cache(int, ST).p̄ .- timestep(int) .* cache(int, ST).f̃ ./ 2
end

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
issymmetric(method::PMVImidpoint) = missing
issymplectic(method::PMVImidpoint) = true


const PMVImidpointIntegrator{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:PMVImidpoint}


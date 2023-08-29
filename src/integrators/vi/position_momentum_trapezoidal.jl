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
issymmetric(method::PMVItrapezoidal) = missing
issymplectic(method::PMVItrapezoidal) = true


const PMVItrapezoidalIntegrator{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:PMVItrapezoidal}


@doc raw"""
Hamilton-Pontryagin Integrator using trapezoidal quadrature.

We consider a discrete Lagrangian of the form
```math
L_d (q_{n}, q_{n+1}) = \frac{h}{2} \big[ L (q_{n}, v_{n+1/2}) + L (q_{n+1}, v_{n+1/2}) \big] ,
```
where $q_{n}$ approximates the solution $q(t_n)$ and $v_{n+1/2}$ is the velocity, which is assumed to be constant in the interval $[t_{n}, t_{n+1}]$.
The discrete Hamilton-Pontryagin action reads
```math
A_d [q_d] = \sum \limits_{n=0}^{N-1} \bigg[ L_d (q_{n}, q_{n+1}) + h \left< p_{n+1} , v_{n+1/2} - \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \right> \bigg] ,
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
0 &= \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n}, v_{n+1/2})
   + \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n+1}, v_{n+1/2}) \\
  &- h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1}
   - h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1} , \\
p_{n+1} &= \frac{1}{2} \, \frac{\partial L}{\partial v} (q_{n}, v_{n+1/2})
         + \frac{1}{2} \, \frac{\partial L}{\partial v} (q_{n+1}, v_{n+1/2}) , \\
v_{n+1/2} &= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) .
\end{aligned}
```
The current implementation requires that $D_2 \phi_h $ and $D_a \phi_h$ do not depend on the second and third argument, but only $D_1 \phi_h$ is allowed to depend on all three arguments.
Otherwise it would be necessary to prescribe initial conditions $(q_{-1}, a_{-1,0}, q_{0}, p_{0})$ instead of just $(q_{0}, p_{0})$.

"""
struct HPItrapezoidal{ϕT, D₁ϕT, D₂ϕT, DₐϕT, PT} <: HPIMethod
    ϕ::ϕT
    D₁ϕ::D₁ϕT
    D₂ϕ::D₂ϕT
    Dₐϕ::DₐϕT
    params::PT
end

isexplicit(method::HPItrapezoidal) = false
isimplicit(method::HPItrapezoidal) = true
issymmetric(method::HPItrapezoidal) = missing
issymplectic(method::HPItrapezoidal) = true


const HPItrapezoidalIntegrator{DT,TT} = GeometricIntegrator{<:Union{IODEProblem{DT,TT},LODEProblem{DT,TT}}, <:HPItrapezoidal}

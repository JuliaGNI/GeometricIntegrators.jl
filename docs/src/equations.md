# Equations

In *GeometricIntegrators.jl* we support three basic types of equations:
* ordinary differential equations (ODEs),
* differential algebraic equations (DAEs),
* stochastic differential equations (SDEs).

For each type, there are several subtypes
* standard equations ([`ODE`](@ref), [`DAE`](@ref), [`SDE`](@ref)),
* implicit equations ([`IODE`](@ref), [`IDAE`](@ref)),
* partitioned equations ([`PODE`](@ref), [`PDAE`](@ref), [`PSDE`](@ref)),
* Hamiltonian equations ([`HODE`](@ref), [`HDAE`](@ref)),
* Lagrangian equations ([`LODE`](@ref), [`LDAE`](@ref)),
* split equations ([`SODE`](@ref), [`SPDAE`](@ref), [`SPSDE`](@ref)).

Each equation holds a number of functions determining the vector field, constraints, initial conditions, and possibly additional information like parameters, periodicity, invariants and the Hamiltonian or Lagrangian.

## Ordinary Differential Equations (ODEs)

Ordinary differential equations define an initial value problem of the form
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``.

### Partitioned Ordinary Differential Equations (PODEs)

A partitioned ordinary differential equation has the form
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Hamiltonian Ordinary Differential Equations (HODEs)

A special case of a partitioned ordinary differential equation is a canonical
Hamiltonian system of equations,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{aligned}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Implicit Ordinary Differential Equations (IODEs)

An implicit ordinary differential equations is of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t))
\end{aligned}
```
with force field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Lagrangian Ordinary Differential Equations (LODEs)

A special case of an implicit ordinary differential equations is a Lagrangian system of equations,
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t))
\end{aligned}
```
with momentum ``p`` and force field ``f``, given by
```math
\begin{aligned}
p &= \frac{\partial L}{\partial v} , &
f &= \frac{\partial L}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the solution ``(q,p)`` taking values
in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.



## Differential Algebraic Equation (DAE)

Differential algebraic initial value problems are of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
0 &= \phi (t, q(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``,
initial conditions ``q_{0}`` and ``\lambda_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}`` and the algebraic variable ``\lambda``
taking values in ``\mathbb{R}^{m}``.


### Partitioned Differential Algebraic Equation (PDAE)

A partitioned differential algebraic equation has the form
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{m}``.

### Hamiltonian Differential Algebraic Equation (HDAE)

A special form of a partitioned differential algebraic equation is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{g}(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{f}(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{m} \times \mathbb{R}^{m}``.


### Implicit Differential Algebraic Equation (IDAE)

An implicit differential algebraic initial value problem takes the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= p(t, q(t), v(t)) , && \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with force field ``f``, the momentum defined by ``p``, projection ``u`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{m}``.


### Lagrangian Differential Algebraic Equation (LDAE)

A special case of an implicit initial value problem is a Lagrangian differential
algebraic equation of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) + \lambda(t), &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), \lambda(t)) + \bar{g} (t, q(t), \mu(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with momentum ``p`` and force field ``f``, given by
```math
\begin{aligned}
p &= \frac{\partial L}{\partial v} , &
f &= \frac{\partial L}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(v, \lambda, \mu)`` taking values in
``\mathbb{R}^{d} \times \mathbb{R}^{m} \times \mathbb{R}^{m}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variables ``v``, ``\lambda`` and ``\mu``.


## Stochastic Differential Equations (SDEs)


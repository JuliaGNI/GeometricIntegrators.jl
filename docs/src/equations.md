# Equation and Problem Types

In *GeometricIntegrators.jl* we support three basic types of equations:

* ordinary differential equations (ODEs),
* differential algebraic equations (DAEs),
* stochastic differential equations (SDEs).

For each type, there are several subtypes:

* standard equations ([`ODEProblem`](@ref), [`DAEProblem`](@ref), [`SDEProblem`](@ref)),
* implicit equations ([`IODEProblem`](@ref), [`IDAEProblem`](@ref)),
* partitioned equations ([`PODEProblem`](@ref), [`PDAEProblem`](@ref), [`PSDEProblem`](@ref)),
* Hamiltonian equations ([`HODEProblem`](@ref), [`HDAEProblem`](@ref)),
* Lagrangian equations ([`LODEProblem`](@ref), [`LDAEProblem`](@ref)),
* split equations ([`SODEProblem`](@ref), [`SPDAEProblem`](@ref), [`SPSDEProblem`](@ref)).

Each equation holds a number of functions determining the vector field, constraints, and possibly additional information like parameters, periodicity, invariants and the Hamiltonian or Lagrangian.

To each equation type there exists a corresponding problem type, which holds the equation, initial conditions, parameters, a time span to integrate over, as well as a time step (which is typically fixed in GeometricIntegrators). In addition, these problem types provide some convenience constructors to generate the equation and the problem at once.

All of the equation and problem types are defined in the [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl) package.
The [GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl) package implements a number of example problems for testing and benchmarking.

## Keyword Arguments

All equation and problem types take the following keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

If not set to their corresponding Null types, the user needs to pass a `NamedTuple` whose values are

* functions for invariants,
* arbitrary data structures for parameters, 
* the same data structure as the solution for periodicity.

The latter should be zero everywhere, except for those components, that are periodic, i.e., whose value are supposed to stay within a range `(0, max)`. Support for ranges starting with other values than zero is currently missing but can be added if demand arises.


## Ordinary Differential Equations (ODEs)

Ordinary differential equations define an initial value problem of the form
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``.

The user needs to specify a function `v` that computes the vector field and must have the interface
```julia
    function v(v, t, q, params)
        v[1] = ...
        v[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
vector which holds the result of evaluating the vector field ``v`` on `t` and
`q`, and `params` are constant parameters on which the vector field may depend.

To create and `ODE`, one only needs to pass this function:
```julia
equ = ODE(v)
```
The full constructor would look like
```julia
equ = ODE(v; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```
where all keyword arguments, namely invariants, parameters and periodicity are by default initialized to be absent.

With this, we can create an `ODEProblem` via
```julia
tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5, 0.0]

prob = GeometricProblem(equ, tspan, tstep, ics = (q=q₀,))
```

Typically, one would create an ODEProblem right away.
We will see this in the next example for the harmonic oscillator.
The dynamical equations are given by
```math
\dot{q} (t) = \begin{pmatrix}
0 & 1 \\
-k & 0 \\
\end{pmatrix} q(t) ,
\qquad
q \in \mathbb{R}^{2} .
```

In order to create an `ODEProblem` for the harmonic oscillator, we need to write the following code:
```julia
function v(v, t, x, params)
    v[1] = x[2]
    v[2] = - params.k * x[1]
end

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5, 0.0]

prob = ODEProblem(v, tspan, tstep, q₀; parameters = (k = 0.5,))
```

The energy of the harmonic oscillator is preserved, so we can add it as an invariant, 
```julia
energy(t, q, params) = q[2]^2 / 2 + params.k * q[1]^2 / 2

prob = ODEProblem(v, tspan, tstep, q₀; parameters = (k = 0.5,), invariants = (h=energy,))
```


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
taking values in ``T^{*} Q \simeq \mathbb{R}^{d} \times \mathbb{R}^{d}``.

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
in ``T^{*} Q \simeq \mathbb{R}^{d} \times \mathbb{R}^{d}``.



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


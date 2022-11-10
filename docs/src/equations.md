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

The latter should be zero everywhere, except for those components, that are periodic, i.e.,
whose value are supposed to stay within a range `(0, max)`. Support for ranges starting
with other values than zero is currently missing but can be added if demand arises.


## Ordinary Differential Equations (ODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.ode_equations)
```

#### Example: Harmonic Oscillator

As an example, let us consider the harmonic oscillator.
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
```@example
using GeometricIntegrators # hide
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

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.pode_equations)
```

#### Example: Harmonic Oscillator

As an example, let us consider the harmonic oscillator.
The dynamical equations are given by
```math
\begin{aligned}
\dot{q} (t) &= p (t) \\
\dot{p} (t) &= - k \, q(t) .
\end{aligned}
```

In order to create a `PODEProblem` for the harmonic oscillator, we need to write the following code:
```@example
using GeometricIntegrators # hide
function v(v, t, q, p, params)
    v[1] = p[1]
end

function f(f, t, q, p, params)
    f[1] = - params.k * q[1]
end

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5]
p₀ = [0.0]

prob = PODEProblem(v, f, tspan, tstep, q₀, p₀; parameters = (k = 0.5,))
```

The energy of the harmonic oscillator is preserved, so we can add it as an invariant, 
```julia
energy(t, q, p, params) = p[1]^2 / 2 + params.k * q[1]^2 / 2

prob = PODEProblem(v, f, tspan, tstep, q₀, p₀; parameters = (k = 0.5,), invariants = (h=energy,))
```


### Hamiltonian Ordinary Differential Equations (HODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.hode_equations)
```

#### Example: Harmonic Oscillator

As an example, let us consider the harmonic oscillator.
The dynamical equations are given by
```math
\begin{aligned}
\dot{q} (t) &= p (t) \\
\dot{p} (t) &= - k \, q(t) ,
\end{aligned}
```
which can also be written as 
```math
\begin{pmatrix}
\dot{q} (t) \\
\dot{p} (t) \\
\end{pmatrix} = \begin{pmatrix}
0 & 1 \\
-1 & 0 \\
\end{pmatrix}
\nabla H( q(t) , p(t) ) ,
\qquad
H(q,p) = \frac{p^2}{2} + k \, \frac{q^2}{2} ,
```
where $H$ is the Hamiltonian, i.e., the total energy of the system.

In order to create a `HODEProblem` for the harmonic oscillator, we need to write the following code:
```@example
using GeometricIntegrators # hide
function v(v, t, q, p, params)
    v[1] = p[1]
end

function f(f, t, q, p, params)
    f[1] = - params.k * q[1]
end

h(t, q, p, params) = p[1]^2 / 2 + params.k * q[1]^2 / 2

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5]
p₀ = [0.0]

prob = HODEProblem(v, f, h, tspan, tstep, q₀, p₀; parameters = (k = 0.5,))
```


### Implicit Ordinary Differential Equations (IODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.iode_equations)
```

#### Example: Harmonic Oscillator

As an example, let us consider the harmonic oscillator.
In implicit form, its equations are given as follows,
```math
\begin{aligned}
\dot{q} (t) &= v(t) , \\
p(t) &= v(t) , \\
\dot{p} (t) &= - k q(t) .
\end{aligned}
```
Here, `v` acts as a Lagrange multiplier that enforces the "constraint" ``p(t) = v(t)``.

In order to create an `IODEProblem` for the harmonic oscillator, we thus need to write the following code:
```@example
using GeometricIntegrators # hide
function p(p, t, q, v, params)
    p[1] = v[1]
end

function f(f, t, q, v, params)
    p[1] = - params.k * q[1]
end

function g(f, t, q, v, params)
    p[1] = - params.k * q[1]
end

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5]
p₀ = [0.0]

prob = IODEProblem(p, f, tspan, tstep, q₀, p₀; parameters = (k = 0.5,))
```

The energy of the harmonic oscillator is preserved, so we can add it as an invariant, 
```julia
energy(t, q, v, params) = v[1]^2 / 2 + params.k * q[1]^2 / 2

prob = IODEProblem(p, f, tspan, tstep, q₀, p₀; parameters = (k = 0.5,), invariants = (h=energy,))
```


### Lagrangian Ordinary Differential Equations (LODEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.lode_equations)
```

#### Example: Harmonic Oscillator

As an example, let us consider the harmonic oscillator.
Its Lagrangian is given by
```math
L(q, \dot{q}) = \frac{\dot{q}^2}{2} - k \, \frac{q^2}{2} ,
```
so that the Euler-Lagrange equations
```math
\frac{d}{dt} \frac{\partial L}{\partial \dot{q}} = \frac{\partial L}{\partial q} ,
```
become
```math
\ddot{q} (t) = - k q(t) .
```

Most integrators for Lagrangian systems do not solve this second order system (semi-spray form),
but instead use a reformulation as an implicit ordinary differential equation.
This formulation can most easily be obtained from a Hamilton-Pontryagin principle
```math
\delta \int \limits_{t_0}^{t_1} \big[ L(q, v) + \left< p , \dot{q} - v \right> \big] = 0 ,
```
as follows,
```math
\begin{aligned}
\dot{q} (t) &= v(t) , \\
p(t) &= \frac{\partial L}{\partial v} (q(t),v(t)) = v(t) , \\
\dot{p} (t) &= \frac{\partial L}{\partial q} (q(t),v(t)) = - k q(t) .
\end{aligned}
```
Here, `v` acts as a Lagrange multiplier that enforces the "constraint" ``p(t) = \partial L / \partial v``.

In order to create an `LODEProblem` for the harmonic oscillator, we thus need to write the following code:
```@example
using GeometricIntegrators # hide
function p(p, t, q, v, params)
    p[1] = v[1]
end

function f(f, t, q, v, params)
    p[1] = - params.k * q[1]
end

function ω(f, t, q, v, params)
    ω[1,1] =  0
    ω[1,2] = -1
    ω[2,1] = +1
    ω[2,2] =  0
end

function l(t, q, v, params)
    v[1]^2 / 2 - params.k * q[1]^2 / 2
end

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5]
p₀ = [0.0]

prob = LODEProblem(p, f, ω, l, tspan, tstep, q₀, p₀; parameters = (k = 0.5,))
```

The energy of the harmonic oscillator is preserved, so we can add it as an invariant, 
```julia
energy(t, q, v, params) = v[1]^2 / 2 + params.k * q[1]^2 / 2

prob = LODEProblem(p, f, ω, l, tspan, tstep, q₀, p₀; parameters = (k = 0.5,), invariants = (h=energy,))
```


## Differential Algebraic Equation (DAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.dae_equations)
```

### Example: Harmonic Oscillator

As an example we consider the harmonic oscillator, with an additional
constraint that enforces energy conservation. While the system itself
is energy conserving, most integrators do not respect this property.
A possible way of remedying this flaw is to explicitly add energy
conservation as an algebraic constraint. 

```@example
using GeometricIntegrators # hide
hamiltonian(t, q, params) = q[2]^2 / 2 + params.k * q[1]^2 / 2

function v(v, t, q, params)
    v[1] = q[2]
    v[2] = - params.k * q[1]
end

function u(u, t, q, λ, params)
    u[1] = λ[1] * params.k * q[1]
    u[2] = λ[1] * q[2]
end

function ϕ(ϕ, t, q, params)
    ϕ[1] = hamiltonian(t, q, params)
end

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [1., 1.]
λ₀ = [0.]
params = (k=0.5,)

dae = DAEProblem(v, u, ϕ, tspan, tstep, q₀, λ₀; parameters = params)
```


### Partitioned Differential Algebraic Equation (PDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.pdae_equations)
```


### Hamiltonian Differential Algebraic Equation (HDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.hdae_equations)
```


### Implicit Differential Algebraic Equation (IDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.idae_equations)
```


### Lagrangian Differential Algebraic Equation (LDAE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.ldae_equations)
```


## Stochastic Differential Equations (SDEs)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.sde_equations)
```

#### Example: Kubo Oscillator

```@example
using GeometricIntegrators # hide
function v(t, q, v, params)
    v[1] = + params.λ * q[2]
    v[2] = - params.λ * q[1]
end

function B(t, q, B, params)
    for j in axes(B, 2)
        B[1,j] = + params.ν * q[2]
        B[2,j] = - params.ν * q[1]
    end
end

tspan = (0.0, 1.0)
tstep = 0.01
q₀ = [0.5, 0.0]

prob = SDEProblem(v, B, tspan, tstep, q₀; parameters = (λ=2.0, μ=1.0))
```


### Partitioned Stochastic Differential Equation (PSDE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.psde_equations)
```

#### Example: Kubo Oscillator

```@example
using GeometricIntegrators # hide
function v(v, t, q, p, params)
    v[1] = + params.λ * p[1]
end

function f(f, t, q, p, params)
    f[1] = - params.λ * q[1]
end

function B(B, t, q, p, params)
    B[1,1] = params.ν * p[1]
end

function G(G, t, q, p, params)
    G[1,1] = - params.ν * q[1]
end

tspan = (0.0, 1.0)
tstep = 0.01
q₀ = [0.5]
p₀ = [0.0]

prob = PSDEProblem(v, f, B, G, tspan, tstep, q₀, p₀; parameters = (λ=2.0, μ=1.0))
```


### Split Partitioned Stochastic Differential Equation (SPSDE)

```@eval
using GeometricEquations, Markdown
Markdown.parse(GeometricEquations.spsde_equations)
```

#### Example: Kubo Oscillator

```@example
using GeometricIntegrators # hide
function v(v, t, q, p, params)
    v[1] = + params.λ * p[1]
end

function f1(f, t, q, p, params)
    f[1] = - params.λ * q[1]
end

function f2(f, t, q, p, params)
    f[1] = 0
end

function B(B, t, q, p, params)
    B[1,1] = params.ν * p[1]
end

function G1(G, t, q, p, params)
    G[1,1] = - params.ν * q[1]
end

function G2(G, t, q, p, params)
    G[1,1] = 0
end

tspan = (0.0, 1.0)
tstep = 0.01
q₀ = [0.5]
p₀ = [0.0]

prob = SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, q₀, p₀; parameters = (λ=2.0, μ=1.0))
```

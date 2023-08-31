```@meta
CurrentModule = GeometricIntegrators
```

# Tutorial

In this tutorial, we try to give an overview of the basic usage of GeometricIntegrators and its main components.


## Installation

*GeometricIntegrators.jl* can be installed using Julia's built-in package manager in the command line interface by
```
julia> ]
(v1.9) pkg> add GeometricIntegrators
```
In a Jupyter notebook, *GeometricIntegrators.jl* can be installed by explicitly using the `Pkg` module as
```julia
using Pkg
Pkg.add("GeometricIntegrators")
```
This will install the library itself as well as all dependencies.


## Basic usage

In the simplest cases, the use of `GeometricIntegrators.jl` requires the
construction of two objects, an equation and an integrator. For many standard
methods, the integrator is implicitly selected by specifying an equation and a
tableau.

Before any use, we need to load GeometricIntegrators,
```@example 1
using GeometricIntegrators
```
Then we can create an ODE object for the equation $\dot{x} (t) = x(t)$ with initial condition $x(0) = 1$, integration time span $(0, 1)$ and a time step of $\Delta t = 0.1$,
```@example 1
prob = ODEProblem((ẋ, t, x, params) -> ẋ[1] = x[1], (0.0, 1.0), 0.1, [1.0])
```
create an integrator for this ODE, using the explicit Euler method
```@example 1
int = GeometricIntegrator(prob, ExplicitEuler())
```
and compute the solution,
```@example 1
sol = integrate(int)
```
Plot and compare with the exact solution
```@example 1
using Plots
plot(xlims=[0,1], xlab="t", ylab="x(t)", legend=:bottomright)
plot!(sol.t, sol.q[:,1], label="numeric")
plot!(sol.t, exp.(sol.t), label="exact")
savefig("images/tutorial-ode-1.png") # hide
```

![](images/tutorial-ode-1.png)


## Equations

In *GeometricIntegrators.jl* we distinguish between three basic types of equations:
* ordinary differential equations (ODEs),
* differential algebraic equations (DAEs),
* stochastic differential equations (SDEs).

For each type, there are several subtypes
* standard equations ([`ODEProblem`](@ref), [`DAEProblem`](@ref), [`SDEProblem`](@ref)),
* implicit equations ([`IODEProblem`](@ref), [`IDAEProblem`](@ref)),
* partitioned equations ([`PODEProblem`](@ref), [`PDAEProblem`](@ref), [`PSDEProblem`](@ref)),
* Hamiltonian equations ([`HODEProblem`](@ref), [`HDAEProblem`](@ref)),
* Lagrangian equations ([`LODEProblem`](@ref), [`LDAEProblem`](@ref)),
* split equations ([`SODEProblem`](@ref), [`SPDAEProblem`](@ref), [`SPSDEProblem`](@ref)).


#### Ordinary differential equations

Consider an ODE of the form
```math
\dot{x} (t) = v(t, x(t)) ,
```
where $\dot{x}$ denotes the derivative of $x$ and $f$ the vector field of the
equation, which is assumed to depend on both $t$ and $x$.
In the following, we will solve the mathematical pendulum, whose equations are
given by
```math
\begin{pmatrix}
\dot{x}_1 \\
\dot{x}_2 \\
\end{pmatrix}
=
\begin{pmatrix}
x_2 \\
\sin (x_1) \\
\end{pmatrix} .
```

Together with the integration time span `(t₀,t₁)` and the time step, an
ODE defines an `ODEProblem`.

The user needs to specify a function `ẋ` that computes the vector field and
must have the interface
```julia
function ẋ(v, t, x, params)
    v[1] = ...
    v[2] = ...
    ...
end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
vector which holds the result of evaluating the vector field ``v`` on `t` and
`q`, and `params` is a `NamedTuple` of constant parameters on which the vector
field may depend.

For the mathematical pendulum, this could look as follows:
```@example 1
function ẋ(v, t, x, params)
    v[1] = x[2]
    v[2] = sin(x[1])
end
```

An `ODEProblem` is instantiated by
```
ODEProblem(<vector field>, <time span>, <time step>, <initial conditions>; kwargs...)
```
so to create and `ODEProblem`, one only needs to pass the above function `ẋ`, a tuple
`tspan` containing the start and end times of the integration, the time step
`tstep` as well as an initial condition:
```@example 1
tspan = (0.0, 10.0)
tstep = 0.1
x₀ = [acos(0.4), 0.0]

ode = ODEProblem(ẋ, tspan, tstep, x₀)
```
The full constructor would look like
```@example 1
ode = ODEProblem(ẋ, tspan, tstep, x₀; invariants = NullInvariants(),
                 parameters = NullParameters(), periodicity = NullPeriodicity())
```
where all keyword arguments, namely invariants, parameters and periodicity, are
by default initialized to be absent.


#### Partitioned ordinary differential equations

The pendulum problem is a Hamiltonian system that can also be expressed as
```math
\dot{q} = \frac{\partial H}{\partial p} = p ,
\hspace{3em}
\dot{p} = - \frac{\partial H}{\partial q} = \sin (q) ,
\hspace{3em}
H (q,p) = \frac{1}{2} p^2 + \cos (q) .
```
This structure, namely the partitioning into two sets of variables $(q,p)$
instead of $x$, can be exploited for more efficient integration.
Such equations can be defined in terms of a partitioned ODE, where the vector
fields are specified separately,
```@example 1
function q̇(v, t, q, p, params)
    v[1] = p[1]
end

function ṗ(f, t, q, p, params)
    f[1] = sin(q[1])
end

pode = PODEProblem(q̇, ṗ, (0.0, 25.0), 0.1, [acos(0.4)], [0.0])
```
The first two arguments to the PODE constructor are the functions that determine
the vector fields of the equations $\dot{q} (t) = v(t, q(t), p(t))$ and
$\dot{p} (t) = f(t, q(t), p(t))$. The third and fourth argument determines the
initial conditions of $q$ and $p$, respectively.
The functions defining the vector field have to take four arguments, the current
time `t`, the current solution vectors `q` and `p` and the output vector
`v` or `f`.


## Integrators

We support a number of standard integrators (geometric and non-geometric) like
explicit, implicit and partitioned Runge-Kutta methods, splitting methods and
general linear methods (_planned_).

In order to instantiate many of the standard integrators, one needs to specify
an ODEProblem, a method and a timestep, e.g.,
```@example 1
int = GeometricIntegrator(ode, ExplicitEuler())
```
In order to run the integrator, the `integrate()` functions is called, passing
an integrator object and the number of time steps to integrate:
```@example 1
sol = integrate(int)
```
The integrate function automatically creates an appropriate solution object,
that contains the result of the integration.

```@example 1
plot(sol.q[:,1], sol.q[:,2], xlab="x(t)", ylab="y(t)", legend=:none)
savefig("images/tutorial-ode-2.png"); nothing # hide
```

![](images/tutorial-ode-2.png)

Observe that the explicit Euler method is not well suited for integrating this
system. The solutions drifts away although it should follow closed orbits.

For a Hamiltonian system, defined as a PODE, a different methods might be more
appropriate, for example a symplectic Euler method,
```@example 1
sol = integrate(pode, LobattoIIIAIIIB(2))
```
This creates a different integrator, which exploits the partitioned structure
of the system. The solution return by the integrate step will also be a different
solution, adapted to the partitioned system.

```@example 1
plot(sol.q[:,1], sol.p[:,1], xlab="q(t)", ylab="p(t)", legend=:none)
savefig("images/tutorial-pode-1.png"); nothing # hide
```

![](images/tutorial-pode-1.png)

Moreover, this method respects the Hamiltonian structure of the system, resulting
in closed orbits following the contours of the system's energy.


## Overview of Available Methods

GeometricIntegrators.jl provides a plethora of geometric integrators as well as non-geometric integrators (mainly for testing and benchmarking purposes).
Most integrators can be selected by a simple method type, which also stores parameters.
Some integrator families can also be selected by specifying a tableau, that is a Butcher tableau for Runge-Kutta methods, a pair of tableaus for partitioned Runge-Kutta and VPRK methods, or generalizations thereof for SPARK methods.
Other integrators, such as Galerkin variational integrators require the specification of a basis and a quadrature rule.

The correct integrator is automatically selected based on the method and problem types by calling
```
GeometricIntegrator(problem, method)
```

As an example, consider an ODE like the harmonic oscillator, which is included in GeometricEquations.jl:
```@example 1
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
```
```@example 1
prob = HarmonicOscillator.odeproblem()
```
Create an explicit Euler method:
```@example 1
method = ExplicitEuler()
```
And now create an Integrator with the general `Integrator` constructor:
```@example 1
int = GeometricIntegrator(prob, method)
```
We see that we obtained an `IntegratorERK`, i.e., an explicit Runge-Kutta integrator.
If instead we choose the implicit Euler method:
```@example 1
method = ImplicitEuler()
```
the general `Integrator` constructor creates a different integrator:
```@example 1
int = GeometricIntegrator(prob, method)
```
namely an `IntegratorFIRK`, i.e., a fully implicit Runge-Kutta integrator.

GeometricIntegrators automatically detects if a Runge-Kutta tableau is explicit, diagonally implicit or fully implicity and creates the corresponding Integrator.

Certain Runge-Kutta method such as Gauß, Radau and Lobatto methods are available for an arbitrary number of stages.
Here the number of stages has to be speficied
```@example 1
int = GeometricIntegrator(prob, Gauss(1))
```
Special integrators, such as Vartiational Partitioned Runge-Kutta (VPRK) methods, can be initialised by providing one or two tableaus, that is
```@example 1
method = VPRK(Gauss(1))
```
or
```@example 1
method = VPRK(LobattoIIIA(2), LobattoIIIB(2))
```
For standard tableaus there also exist shortcuts, such as
```@example 1
method = VPRKGauss(1)
```
or
```@example 1
method = VPRKLobattoIIIAIIIB(2)
```
For the purpose of a complete example, consider again the harmonic oscillator:
```@setup 2
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
```
```example 2
prob = HarmonicOscillator.iodeproblem()
```
Create a VPRK tableau that uses Gauss-Legendre Runge-Kutta coefficients with two stages:
```example 2
method = VPRKGauss(2)
```
If we call the `Integrator` constructor,
```example 2
int = Integrator(prob, method)
```
we obtain a `IntegratorVPRK`.

Once an integrator is obtained, we can just call the function
```
integrate(integrator)
```
to perform the actual integration steps, where `ntime` defines the number of steps to integrate:
```@setup 3
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
prob = HarmonicOscillator.odeproblem()
```
```@example 3
int = GeometricIntegrator(prob, ExplicitEuler())
sol = integrate(int)
```
The `integrate` function returns a solution object that stores the solution for each time step.
If the solution object is created manually, there exists a function
```
integrate!(integrator, solution)
```
that operates on an existing solution.


### Integrators for ODEs

The main method types for ODEs currently implemented are Runge-Kutta methods and splitting methods.

#### Runge-Kutta Methods

Any Runge-Kutta method can be selected by the [`RK`](@ref) method
```julia
rk = RK(tableau)
```
where `tableau` is any tableau from [RungeKutta.Tableaus](@ref RungeKuttaTableaus). For most tableaus there also exist explicit shortcuts to select the method. These are listed in the following.


#### Explicit Runge-Kutta Methods

| Function                        | Order | Method                           |
|:--------------------------------|:------|:---------------------------------|
| [`ExplicitEuler`](@ref)         | 1     | Explicit / Forward Euler         |
| [`ExplicitMidpoint`](@ref)      | 2     | Explicit Midpoint                |
| [`Heun2`](@ref)                 | 2     | Heun's Method of order two       |
| [`Heun3`](@ref)                 | 3     | Heun's Method of order three     |
| [`Ralston2`](@ref)              | 2     | Ralston's Method of order two    |
| [`Ralston3`](@ref)              | 3     | Ralston's Method of order three  |
| [`Runge2`](@ref)                | 2     | Runge's Method                   |
| [`Kutta3`](@ref)                | 3     | Kutta's Method                   |
| [`RK416`](@ref)                 | 4     | Explicit 4th order Runge-Kutta (1/6 rule) |
| [`RK438`](@ref)                 | 4     | Explicit 4th order Runge-Kutta (3/8 rule) |


#### Diagonally Implicit Runge-Kutta Methods

| Function                        | Order | Method                           |
|:--------------------------------|:------|:---------------------------------|
| [`CrankNicolson`](@ref)         | 3     | Crank-Nicholson Method           |
| [`KraaijevangerSpijker`](@ref)  | 3     | Kraaijevanger & Spijker's Method |
| [`QinZhang`](@ref)              | 3     | Qin & Zhang's Method             |
| [`Crouzeix`](@ref)              | 3     | Crouzeix's Method                |


#### Fully Implicit Runge-Kutta Methods

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`ImplicitEuler`](@ref)         | 1     | Implicit / Backward Euler   |
| [`ImplicitMidpoint`](@ref)      | 2     | Implicit Midpoint           |
| [`SRK3`](@ref)                  | 4     | Symmetric Runge-Kutta s=3   |


#### Gauß, Radau and Lobatto Methods

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`Gauss`](@ref)                 | 2s    | Gauss-Legendre              |
| [`RadauIA`](@ref)               | 2s-1  | Radau-IA                    |
| [`RadauIB`](@ref)               | 2s-1  | Radau-IB                    |
| [`RadauIIA`](@ref)              | 2s-1  | Radau-IIA                   |
| [`RadauIIB`](@ref)              | 2s-1  | Radau-IIB                   |
| [`LobattoIII`](@ref)            | 2s-2  | Lobatto-III                 |
| [`LobattoIIIA`](@ref)           | 2s-2  | Lobatto-IIIA                |
| [`LobattoIIIB`](@ref)           | 2s-2  | Lobatto-IIIB                |
| [`LobattoIIIC`](@ref)           | 2s-2  | Lobatto-IIIC                |
| [`LobattoIIID`](@ref)           | 2s-2  | Lobatto-IIID                |
| [`LobattoIIIE`](@ref)           | 2s-2  | Lobatto-IIIE                |
| [`LobattoIIIF`](@ref)           | 2s    | Lobatto-IIIF                |
| [`LobattoIIIF`](@ref)           | 2s    | Lobatto-IIIF                |
| [`LobattoIIIG`](@ref)           | 2s    | Lobatto-IIIG                |

All of these tableaus are generated on the fly and take the number of stages `s` as parameter.


#### Splitting Methods

| Function                        | Order | Method                       |
|:--------------------------------|:------|:-----------------------------|
| [`LieA`](@ref)                  | 1     | Lie-Trotter Splitting A      |
| [`LieB`](@ref)                  | 1     | Lie-Trotter Splitting B      |
| [`Strang`](@ref)                | 2     | Strang / Marchuk Splitting   |
| [`Marchuk`](@ref)               | 2     | Strang / Marchuk Splitting   |
| [`StrangA`](@ref)               | 2     | Strang / Marchuk Splitting A |
| [`StrangB`](@ref)               | 2     | Strang / Marchuk Splitting B |
| [`McLachlan2`](@ref)            | 2     | McLachlan's 2nd order symmetric, minimum error composition method |
| [`McLachlan4`](@ref)            | 2     | McLachlan's 4th order symmetric, minimum error composition method |
| [`TripleJump`](@ref)            | 4     | 4th order "Triple Jump" composition method                        |
| [`SuzukiFractal`](@ref)         | 4     | Suzuki's 4th order "fractal" composition method                   |


### Integrators for partitioned ODEs

#### Partitioned Runge-Kutta Methods

Any partitioned Runge-Kutta method can be selected by the [`PRK`](@ref) method
```julia
prk = PRK(tableau)
```
where `tableau` is any tableau from [RungeKutta.PartitionedTableaus](@ref partitioned-runge-kutta-rableaus).
For most tableaus there also exist explicit shortcuts to select the method. These are listed in the following.

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`LobattoIIIAIIIB`](@ref)       | 2s-2  | Lobatto-IIIA-IIIB           |
| [`LobattoIIIBIIIA`](@ref)       | 2s-2  | Lobatto-IIIB-IIIA           |
| [`LobattoIIIAIIIĀ`](@ref)       | 2s-2  | Lobatto-IIIA-IIIĀ           |
| [`LobattoIIIBIIIB̄`](@ref)       | 2s-2  | Lobatto-IIIB-IIIB̄           |
| [`LobattoIIICIIIC̄`](@ref)       | 2s-2  | Lobatto-IIIC-IIIC̄           |
| [`LobattoIIIC̄IIIC`](@ref)       | 2s-2  | Lobatto-IIIC̄-IIIC           |
| [`LobattoIIIDIIID̄`](@ref)       | 2s-2  | Lobatto-IIID-IIID̄           |
| [`LobattoIIIEIIIĒ`](@ref)       | 2s-2  | Lobatto-IIIE-IIIĒ           |
| [`LobattoIIIFIIIF̄`](@ref)       | 2s    | Lobatto-IIIF-IIIF̄           |
| [`LobattoIIIF̄IIIF`](@ref)       | 2s    | Lobatto-IIIF̄-IIIF           |
| [`LobattoIIIGIIIḠ`](@ref)       | 2s    | Lobatto-IIIG-IIIḠ           |


### Integrators for implicit ODEs

All implicit Runge-Kutta and partitioned Runge-Kutta methods can also be applied to implicit ODEs.

| Function                        | Order | Method                      |
|:--------------------------------|:------|:----------------------------|
| [`ImplicitEuler`](@ref)         | 1     | Implicit / Backward Euler   |
| [`ImplicitMidpoint`](@ref)      | 2     | Implicit Midpoint           |
| [`SRK3`](@ref)                  | 4     | Symmetric Runge-Kutta s=3   |
|                                 |       |                             |
| [`Gauss`](@ref)                 | 2s    | Gauss-Legendre              |
| [`RadauIA`](@ref)               | 2s-1  | Radau-IA                    |
| [`RadauIB`](@ref)               | 2s-1  | Radau-IB                    |
| [`RadauIIA`](@ref)              | 2s-1  | Radau-IIA                   |
| [`RadauIIB`](@ref)              | 2s-1  | Radau-IIB                   |
| [`LobattoIII`](@ref)            | 2s-2  | Lobatto-III                 |
| [`LobattoIIIA`](@ref)           | 2s-2  | Lobatto-IIIA                |
| [`LobattoIIIB`](@ref)           | 2s-2  | Lobatto-IIIB                |
| [`LobattoIIIC`](@ref)           | 2s-2  | Lobatto-IIIC                |
| [`LobattoIIID`](@ref)           | 2s-2  | Lobatto-IIID                |
| [`LobattoIIIE`](@ref)           | 2s-2  | Lobatto-IIIE                |
| [`LobattoIIIF`](@ref)           | 2s    | Lobatto-IIIF                |
| [`LobattoIIIG`](@ref)           | 2s    | Lobatto-IIIG                |
|                                 |       |                             |
| [`LobattoIIIAIIIB`](@ref)       | 2s-2  | Lobatto-IIIA-IIIB           |
| [`LobattoIIIBIIIA`](@ref)       | 2s-2  | Lobatto-IIIB-IIIA           |
| [`LobattoIIIAIIIĀ`](@ref)       | 2s-2  | Lobatto-IIIA-IIIĀ           |
| [`LobattoIIIBIIIB̄`](@ref)       | 2s-2  | Lobatto-IIIB-IIIB̄           |
| [`LobattoIIICIIIC̄`](@ref)       | 2s-2  | Lobatto-IIIC-IIIC̄           |
| [`LobattoIIIC̄IIIC`](@ref)       | 2s-2  | Lobatto-IIIC̄-IIIC           |
| [`LobattoIIIDIIID̄`](@ref)       | 2s-2  | Lobatto-IIID-IIID̄           |
| [`LobattoIIIEIIIĒ`](@ref)       | 2s-2  | Lobatto-IIIE-IIIĒ           |
| [`LobattoIIIFIIIF̄`](@ref)       | 2s    | Lobatto-IIIF-IIIF̄           |
| [`LobattoIIIF̄IIIF`](@ref)       | 2s    | Lobatto-IIIF̄-IIIF           |
| [`LobattoIIIGIIIḠ`](@ref)       | 2s    | Lobatto-IIIG-IIIḠ           |


### Integrators for Lagrangian ODEs

Regular (non-degenerate) Lagragian ODEs can be integrated with Variational Partitioned Runge-Kutta ([`VPRK`](@ref))
methods or Continuous Galerkin Variational Integrators ([`CGVI`](@ref)).

| Function                        | Method                                                             |
|:--------------------------------|:-------------------------------------------------------------------|
| [`VPRK`](@ref)                  | Variational Partitioned Runge-Kutta integrator                     |
|                                 |                                                                    |
| [`VPRKGauss`](@ref)             | VPRK integrator with [`Gauss`](@ref)                               |
| [`VPRKRadauIIA`](@ref)          | VPRK integrator with [`RadauIIA`](@ref)                            |
| [`VPRKRadauIIB`](@ref)          | VPRK integrator with [`RadauIIB`](@ref)                            |
| [`VPRKLobattoIII`](@ref)        | VPRK integrator with [`LobattoIII`](@ref)                          |
| [`VPRKLobattoIIIA`](@ref)       | VPRK integrator with [`LobattoIIIA`](@ref)                         |
| [`VPRKLobattoIIIB`](@ref)       | VPRK integrator with [`LobattoIIIB`](@ref)                         |
| [`VPRKLobattoIIIC`](@ref)       | VPRK integrator with [`LobattoIIIC`](@ref)                         |
| [`VPRKLobattoIIID`](@ref)       | VPRK integrator with [`LobattoIIID`](@ref)                         |
| [`VPRKLobattoIIIE`](@ref)       | VPRK integrator with [`LobattoIIIE`](@ref)                         |
| [`VPRKLobattoIIIF`](@ref)       | VPRK integrator with [`LobattoIIIF`](@ref)                         |
| [`VPRKLobattoIIIG`](@ref)       | VPRK integrator with [`LobattoIIIG`](@ref)                         |
| [`VPRKLobattoIIIAIIIB`](@ref)   | VPRK integrator with [`LobattoIIIAIIIB`](@ref)                     |
| [`VPRKLobattoIIIBIIIA`](@ref)   | VPRK integrator with [`LobattoIIIBIIIA`](@ref)                     |
| [`VPRKLobattoIIIAIIIĀ`](@ref)   | VPRK integrator with [`LobattoIIIAIIIĀ`](@ref)                     |
| [`VPRKLobattoIIIBIIIB̄`](@ref)   | VPRK integrator with [`LobattoIIIBIIIB̄`](@ref)                     |
| [`VPRKLobattoIIICIIIC̄`](@ref)   | VPRK integrator with [`LobattoIIICIIIC̄`](@ref)                     |
| [`VPRKLobattoIIIC̄IIIC`](@ref)   | VPRK integrator with [`LobattoIIIC̄IIIC`](@ref)                     |
| [`VPRKLobattoIIIDIIID̄`](@ref)   | VPRK integrator with [`LobattoIIIDIIID̄`](@ref)                     |
| [`VPRKLobattoIIIEIIIĒ`](@ref)   | VPRK integrator with [`LobattoIIIEIIIĒ`](@ref)                     |
| [`VPRKLobattoIIIFIIIF̄`](@ref)   | VPRK integrator with [`LobattoIIIFIIIF̄`](@ref)                     |
| [`VPRKLobattoIIIF̄IIIF`](@ref)   | VPRK integrator with [`LobattoIIIF̄IIIF`](@ref)                     |
| [`VPRKLobattoIIIGIIIḠ`](@ref)   | VPRK integrator with [`LobattoIIIGIIIḠ`](@ref)                     |


### Integrators for Degenerate Lagrangian ODEs

Degenerate Lagragian ODEs can be integrated with [Degenerate Variational Integrators](integrators/dvi.md) (see also [`DVRK`](@ref))
or Projected Variational Partitioned Runge-Kutta ([`ProjectedVPRK`](@ref)) methods.

| Function                        | Method                                                                    |
|:--------------------------------|:--------------------------------------------------------------------------|
| [`DVIA`](@ref)                  | Symplectic Euler-A Degenerate Variational Integrator                      |
| [`DVIB`](@ref)                  | Symplectic Euler-B Degenerate Variational Integrator                      |
| [`CMDVI`](@ref)                 | Midpoint Degenerate Variational Integrator                                |
| [`CTDVI`](@ref)                 | Trapezoidal Degenerate Variational Integrator                             |
| [`DVRK`](@ref)                  | Degenerate Variational Runge-Kutta integrator                             |
| [`ProjectedVPRK`](@ref)         | Projected Variational Partitioned Runge-Kutta integrator                  |
|                                 |                                                                           |
| [`VPRKpInternal`](@ref)         | VPRK integrator with projection on internal stages                        |
| [`VPRKpLegendre`](@ref)         | VPRK integrator with Legendre projection                                  |
| [`VPRKpMidpoint`](@ref)         | VPRK integrator with Midpoint projection                                  |
| [`VPRKpSecondary`](@ref)        | VPRK integrator with projection on secondary constraint                   |
| [`VPRKpStandard`](@ref)         | VPRK integrator with standard projection                                  |
| [`VPRKpSymmetric`](@ref)        | VPRK integrator with symmetric projection                                 |
| [`VPRKpSymplectic`](@ref)       | VPRK integrator with symplectic projection                                |
| [`VPRKpVariational`](@ref)      | VPRK integrator with variational projection                               |
| [`VPRKpVariationalP`](@ref)     | VPRK integrator with variational projection on P                          |
| [`VPRKpVariationalQ`](@ref)     | VPRK integrator with variational projection on Q                          |


### Integrators for DAEs


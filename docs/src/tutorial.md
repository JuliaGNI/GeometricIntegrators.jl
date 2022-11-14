```@meta
CurrentModule = GeometricIntegrators
```

# Tutorial

In this tutorial, we try to give an overview of the basic usage of GeometricIntegrators and its main components.


## Installation

*GeometricIntegrators.jl* can be installed using Julia's built-in package manager in the command line interface by
```julia
julia> ]
(v1.8) pkg> add GeometricIntegrators
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
create an integrator for this ODE, using the tableau for the explicit Euler method
```@example 1
int = Integrator(prob, TableauExplicitEuler())
```
and compute the solution,
```@example 1
sol = integrate(prob, int)
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
an ODEProblem, a tableau and a timestep, e.g.,
```@example 1
int = Integrator(ode, TableauExplicitEuler())
```
In order to run the integrator, the `integrate()` functions is called, passing
an integrator object and the number of time steps to integrate:
```@example 1
sol = integrate(ode, int)
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

For a Hamiltonian system, defined as a PODE, a different tableau might be more
appropriate, for example a symplectic Euler method,
```@example 1
int = Integrator(pode, TableauLobattoIIIAIIIB(2))
sol = integrate(pode, int)
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


## Tableaus

Many tableaus for Runge-Kutta methods are predefined and can easily be used
like outlined above. For an overview see [here](integrators/overview.md).


#### Custom Tableaus

If required, it is straight-forward to create a custom tableau.
The tableau of Heun's method, for example, is defined as follows:
```@example 1
a = [[0.0 0.0]
     [1.0 0.0]]
b = [0.5, 0.5]
c = [0.0, 1.0]
o = 2

tab = Tableau(:heun, o, a, b, c)
```
Here, `o` is the order of the method, `a` are the coefficients, `b` the weights
and `c` the nodes. For partitioned Runge-Kutta tableaus, `PartitionedTableau` can
be used. The first parameter of the constructor of each tableau assigns a name to
the tableau.
Such custom tableaus can be used in exactly the same as standard tableaus, e.g., by
```@example 1
int = Integrator(ode, tab)
sol = integrate(ode, int)
```
making it very easy to implement and test new methods.


## Solutions

In what we have seen so far, the solution was always automatically created by
the `integrate()` function. While this is often convenient, it is sometimes not
performant, e.g., when carrying out long-time simulations with intermediate
saving of the solution.
In such cases, it is better to preallocate a solution object by
```@example 1
sol = Solution(ode)
```
where the first argument is an equation, the second argument is the time step
and the third argument is the number of time steps that will be computed in one
integration step.
The call to the integrator is then made via
```@example 1
integrate!(int, sol)
```
If several integration cycles shall be performed, the `reset!()` function can be
used to copy the solution of the last time step to the initial conditions of the
solution,
```julia
for i in 1:10
    # integrate!(int, sol)
    #
    # save or process solution
    #
    # reset!(sol)
end
```
All solutions have a `t` field holding the series of time steps that has been
computed in addition to several data fields, for example `q` for an ODE solution,
`q` and `p` for a PODE solution, `q`and `λ` for a DAE solution, and `q`, `p` and
`λ` for a PDAE solution.

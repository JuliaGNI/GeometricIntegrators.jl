
# Tutorial

In the simplest cases, the use of `GeometricIntegrators.jl` requires the
construction of two objects, an equation and an integrator. The integrator
is usually implicitly selected by specifying an equation and a tableau.


## Equations

In *GeometricIntegrators.jl* we distinguish between three basic types of equations:
ordinary differential equations (ODEs), differential algebraic equations (DAEs)
and stochastic differential equations (SDEs).
For each type, there are several subtypes like implicit equations (IODE, etc.),
partitioned equations (PODE, etc.) or split equations (SODE, etc.).

Instantiating an ODE object for the pendulum problem
\\[
\dot{x}_1 = x_2 ,
\hspace{3em}
\dot{x}_2 = \sin (x_1) ,
\\]
can be achieved by
```julia
function pendulum_rhs(t, x, f)
    f[1] = x[2]
    f[2] = sin(x[1])
end

ode = ODE(pendulum_rhs, [acos(0.4), 0.0])
```
The first argument to the ODE constructor is the function that determines the
vector field of the equation $\dot{x} (t) = f(t, x(t))$, and the second argument
determines the initial conditions.
The function defining the vector field has to take three arguments, the current
time `t`, the current solution vector `x` and the output vector `f`.

The pendulum problem is a Hamiltonian system that can also be expressed as
\\[
\dot{q} = \frac{\partial H}{\partial p} = p ,
\hspace{3em}
\dot{p} = - \frac{\partial H}{\partial q} = \sin (q) ,
\hspace{3em}
H (q,p) = \frac{1}{2} p^2 + \cos (q) .
\\]
This structure, namely the partitioning into two sets of variables $(q,p)$
instead of $x$, can be exploited for more efficient integration.
Such equations can be defined in terms of a partitioned ODE, where the vector
fields are specified separately,
```julia
function pendulum_v(t, q, p, v)
    v[1] = p[1]
end

function pendulum_f(t, q, p, f)
    f[1] = sin(q[1])
end

pode = PODE(pendulum_v, pendulum_f, [acos(0.4)], [0.0])
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
an ODE, a tableau and a timestep, e.g.,
```julia
int = Integrator(ode, getTableauExplicitEuler(), 0.1)
```
In order to run the integrator, the `integrate()` functions is called, passing
an integrator object and the number of time steps to integrate:
```julia
sol = integrate(int, 10)
```
The integrate function automatically creates an appropriate solution object,
that contains the result of the integration.

For a Hamiltonian system, defined as a PODE, a different tableau might be more
appropriate, for example a symplectic Euler method,
```julia
int = Integrator(pode, getTableauSymplecticEulerA(), 0.1)
sol = integrate(int, 10)
```
This creates a different integrator, which exploits the partitioned structure
of the system. The solution return by the integrate step will also be a different
solution, adapted to the partitioned system.


## Tableaus

Many tableaus for Runge-Kutta methods are predefined and can easily be used
like outlined above. In particular, this includes the following methods:

#### Explicit Runge-Kutta Methods

| Function           | Order | Method             |
|--------------------|-------|--------------------|
| `getTableauExplicitEuler()`    | 1 | Explicit / Forward Euler    |
| `getTableauExplicitMidpoint()` | 2 | Explicit Midpoint |
| `getTableauHeun()`   | 2 | Heun's Method  |
| `getTableauKutta()`  | 3 | Kutta's Method |
| `getTableauERK4()`   | 4 | Explicit 4th order Runge-Kutta (1/6 rule) |
| `getTableauERK438()` | 4 | Explicit 4th order Runge-Kutta (3/8 rule) |


#### Fully Implicit Runge-Kutta Methods

| Function           | Order | Method             |
|--------------------|-------|--------------------|
| `getTableauImplicitEuler()`    | 1  | Implicit / Backward Euler  |
| `getTableauImplicitMidpoint()` | 2  | Implicit Midpoint          |
| `getTableauRadIIA2()`          | 3  | Radau-IIA s=2              |
| `getTableauRadIIA3()`          | 5  | Radau-IIA s=3              |
| `getTableauSRK3()`             | 4  | Symmetric Runge-Kutta s=3  |
| `getTableauGLRK(s)`            | 2s | Gauss-Legendre Runge-Kutta |


#### Explicit Partitioned Runge-Kutta Methods

| Function           | Order | Method             |
|--------------------|-------|--------------------|
| `getTableauSymplecticEulerA()` | 1 | Symplectic Euler A |
| `getTableauSymplecticEulerB()` | 1 | Symplectic Euler B |
| `getTableauLobattoIIIAIIIB2()` | 2 | Lobatto-IIIA-IIIB  |
| `getTableauLobattoIIIBIIIA2()` | 2 | Lobatto-IIIB-IIIA  |


#### Custom Tableaus

If required, it is straight-forward to create a custom tableau.
The tableau of Heun's method, for example, is defined as follows:
```julia
a = [[0.0 0.0]
     [1.0 0.0]]
b = [0.5, 0.5]
c = [0.0, 1.0]
o = 2

tab = TableauERK(:heun, o, a, b, c)
```
Here, `o` is the order of the method, `a` are the coefficients, `b` the weights
and `c` the nodes. `TableauERK` states that the method is explicit. Other choices
include `TableauFIRK` for fully implicit Runge-Kutta methods, `TableauDIRK` for
diagonally implicit and `TableauSIRK` for singly implicit Runge-Kutta methods.
`TableauEPRK` and `TableauIPRK` can be used for explicit and implicit partitioned
Runge-Kutta methods. The first parameter of the constructor of each tableau
assigns a name to the tableau.
Such custom tableaus can be used in exactly the same as standard tableaus, e.g., by
```julia
int = Integrator(ode, tab, 0.1)
sol = integrate(int, 10)
```
making it very easy to implement and test new methods.


## Solutions

In what we have seen so far, the solution was always automatically created by
the `integrate()` function. While this is often convenient, it is sometimes not
performant, e.g., when carrying out long-time simulations with intermediate
saving of the solution.
In such cases, it is better to preallocate a solution object by
```julia
sol = Solution(ode, 0.1, 10)
```
where the first argument is an equation, the second argument is the time step
and the third argument is the number of time steps that will be computed in one
integration step.
The call to the integrator is then made via
```julia
integrate!(int, sol)
```
If several integration cycles shall be performed, the `reset!()` function can be
used to copy the solution of the last time step to the initial conditions of the
solution,
```julia
for i in 1:10
    integrate!(int, sol)
    #
    # save or process solution
    #
    reset!(sol)
end
```
All solutions have a `t` field holding the series of time steps that has been
computed in addition to several data fields, for example `q` for an ODE solution,
`q` and `p` for a PODE solution, `q`and `λ` for a DAE solution, and `q`, `p` and
`λ` for a PDAE solution.

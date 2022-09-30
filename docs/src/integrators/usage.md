```@meta
CurrentModule = GeometricIntegrators
```

# Integrators

GeometricIntegrators.jl provides a plethora of geometric integrators as well as non-geometric integrators (mainly for testing and benchmarking purposes).
Most integrators are specified by a tableau, that is a Butcher tableau for Runge-Kutta methods, a pair of tableaus for partitioned Runge-Kutta and VPRK methods, or generalizations thereof for SPARK methods.
Other integrators, such as Galerkin variational integrators require the specification of a basis and a quadrature rule.

In many cases, the correct integrator is automatically selected based on the tableau and problem types by calling
```
Integrator(problem, tableau)
```
where `Δt` is the time step.

As an example, consider an ODE like the harmonic oscillator, which is included in GeometricProblems.jl:
```@example 1
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
```
```@example 1
prob = harmonic_oscillator_ode()
```
Create an explicit Euler tableau:
```@example 1
tab = TableauExplicitEuler()
```
And now create an Integrator with the general `Integrator` constructor:
```@example 1
int = Integrator(prob, tab)
```
We see that we obtained an `IntegratorERK`, i.e., an explicit Runge-Kutta integrator.
If instead we choose the implicit Euler tableau:
```@example 1
tab = TableauImplicitEuler()
```
the general `Integrator` constructor creates a different integrator:
```@example 1
int = Integrator(prob, tab)
```
namely an `IntegratorFIRK`, i.e., a fully implicit Runge-Kutta integrator.

This is possible because most integrators come with a dedicated tableau type, so that `Integrator` can dispatch on that.


In some cases, in particular the VPRK integrators, the integrator has to be explicitly specified as there are different integrators that use the same tableau type and operate on the same equation type, here `TableauVPRK` and `IODE`.
Consider again the harmonic oscillator:
```@setup 2
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
```
```example 2
prob = harmonic_oscillator_iode()
```
Create a VPRK tableau that uses Gauss-Legendre Runge-Kutta coefficients with two stages:
```example 2
tab = TableauVPGLRK(2)
```
If we just call the `Integrator` constructor,
```example 2
int = Integrator(prob, tab)
```
we obtain a plain `IntegratorVPRK`.
If we want to use any of the projection methods, we have to explicitly specify the corresponding integrator type:
```example 2
using GeometricIntegrators.Integrators.VPRK
int = IntegratorVPRKpStandard(prob, tab)
```
or
```example 2
int = IntegratorVPRKpSymmetric(prob, tab)
```


Once an integrator is obtained, we can just call the function
```
integrate(problem, integrator)
```
to perform the actual integration steps, where `ntime` defines the number of steps to integrate:
```@setup 3
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
prob = harmonic_oscillator_ode()
```
```@example 3
tab = TableauExplicitEuler()
int = Integrator(prob, tab)
sol = integrate(prob, int)
```
The `integrate` function returns a solution object that stores the solution for each time step.
If the solution object is created manually, there exists a function
```
integrate!(integrator, solution)
```
that operates on an existing solution.

```@meta
CurrentModule = GeometricIntegrators
```

# Integrators

GeometricIntegrators.jl provides a plethora of geometric and non-geometric integrators.
Most integrators are specified by a tableau, that is a Butcher tableau for Runge-Kutta methods, a pair of tableaus for partitioned Runge-Kutta and VPRK methods, or generalizations thereof for SPARK methods.
Other integrators, such as Galerkin variational integrators require the specification of a basis and a quadrature rule.

In many cases, the correct integrator is automatically selected based on the tableau and equation types by calling
```
Integrator(equation, tableau, Δt)
```
where `Δt` is the time step.

As an example, consider an ODE like the harmonic oscillator, which is included as an example problem:
```@setup 1
using GeometricIntegrators
```
```@example 1
ode = TestProblems.HarmonicOscillatorProblem.harmonic_oscillator_ode()
```
Create an explicit Euler tableau:
```@example 1
tab = TableauExplicitEuler()
```
And now create an Integrator with the general `Integrator` constructor:
```@example 1
int = Integrator(ode, tab, 0.1)
```
We see that we obtained an `IntegratorERK`, i.e., an explicit Runge-Kutta integrator.
If instead we choose the implicit Euler tableau:
```@example 1
tab = TableauImplicitEuler()
```
the general `Integrator` constructor creates a different integrator:
```@example 1
int = Integrator(ode, tab, 0.1)
```
namely an `IntegratorFIRK`, i.e., a fully implicit Runge-Kutta integrator.

This is possible because most integrators come with a dedicated tableau type, so that `Integrator` can dispatch on that.


In some cases, in particular the VPRK integrators, the integrator has to be explicitly specified as there are different integrators that use the same tableau type and operate on the same equation type, here `TableauVPRK` and `IODE`.
Consider again the harmonic oscillator:
```@setup 2
using GeometricIntegrators
```
```example 2
iode = harmonic_oscillator_iode
```
Create a VPRK tableau that uses Gauss-Legendre Runge-Kutta coefficients with two stages:
```example 2
tab = TableauVPGLRK(2)
```
If we just call the `Integrator` constructor,
```example 2
int = Integrator(iode, tab, 0.1)
```
we obtain a plain `IntegratorVPRK`.
If we want to use any of the projection methods, we have to explicitly specify the corresponding integrator type:
```example 2
int = IntegratorVPRKpStandard(iode, tab, 0.1)
```
or
```example 2
int = IntegratorVPRKpSymmetric(iode, tab, 0.1)
```


Once an integrator is obtained, we can just call the function
```
integrate(equation, integrator, ntime)
```
to perform the actual integration steps, where `ntime` defines the number of steps to integrate:
```@setup 3
using GeometricIntegrators
ode = TestProblems.HarmonicOscillatorProblem.harmonic_oscillator_ode()
```
```@example 3
tab = TableauExplicitEuler()
int = Integrator(ode, tab, 0.1)
sol = integrate(ode, int, 100)
```
The `integrate` function returns a solution object that stores the solution for each of the `ntime` time steps.
There is also a convenience function that combines all of the above steps in one single call, namely
```
integrate(equation, tableau, Δt, ntime)
```
If the solution object is created manually, there exists a function
```
integrate!(integrator, solution)
```
that operates on an existing solution.

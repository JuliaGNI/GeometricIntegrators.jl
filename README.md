
# GeometricIntegrators.jl

*Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.*

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/GeometricIntegrators.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliagni.github.io/GeometricIntegrators.jl/latest/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
![CI](https://github.com/JuliaGNI/GeometricIntegrators.jl/workflows/CI/badge.svg)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricIntegrators.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricIntegrators.html)
[![codecov Status](https://codecov.io/gh/JuliaGNI/GeometricIntegrators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/GeometricIntegrators.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3648325.svg)](https://doi.org/10.5281/zenodo.3648325)

GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles.  


## Installation

*GeometricIntegrators.jl* is registered in the Julia package registry and can thus easily be installed in the Julia REPL by typing 
```
]add GeometricIntegrators
```


## Basic usage

In the simplest cases, the use of *GeometricIntegrators.jl* requires the
construction of two objects, an equation and an integrator. For many standard
methods, the integrator can be implicitly selected by specifying an equation
and a tableau.

Before any use, we need to load `GeometricIntegrators`,
```julia
using GeometricIntegrators
```
Then we can create an `ODE` object for the equation ẋ(t) = x(t) with initial condition x(0)=1 by
```julia
ode = ODE((t, x, ẋ) -> ẋ[1] = x[1], [1.0]);
```
An integrator for this ODE, using the tableau for the explicit Euler method and a time step of Δt=0.1, is obtained by
```julia
int = Integrator(ode, TableauExplicitEuler(), 0.1);
```
With that, the solution for nₜ=10 time steps is simply computed by
```julia
sol = integrate(ode, int, 10);
```
which returns an object holding the solution for all time steps.
With the help of the *[Plots.jl](https://github.com/JuliaPlots/Plots.jl)* package we can visualise the result and compare with the exact solution:
```julia
using Plots
plot(xlims=[0,1], xlab="t", ylab="x(t)", legend=:bottomright)
plot!(sol.t, sol.q[1,:], label="numeric")
plot!(sol.t, exp.(sol.t), label="exact")
```

![expontential_growth](https://user-images.githubusercontent.com/21168502/100005439-3e30e400-2dc9-11eb-97a9-d485f4e56d86.png)


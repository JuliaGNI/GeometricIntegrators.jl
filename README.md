
# GeometricIntegrators.jl

*Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.*

[![License](https://img.shields.io/badge/license-MIT%20License-blue.svg)](LICENSE.md)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://ddmgni.github.io/GeometricIntegrators.jl/latest/)
[![Build Status](https://travis-ci.org/DDMGNI/GeometricIntegrators.jl.svg?branch=master)](https://travis-ci.org/DDMGNI/GeometricIntegrators.jl)
[![Coverage Status](https://coveralls.io/repos/github/DDMGNI/GeometricIntegrators.jl/badge.svg)](https://coveralls.io/github/DDMGNI/GeometricIntegrators.jl)
[![codecov](https://codecov.io/gh/DDMGNI/GeometricIntegrators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/DDMGNI/GeometricIntegrators.jl)

GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles.  

*Disclaimer:* The package is currently under heavy development. Regular refactoring and breaking changes are expected.


## Features

The following list provides and overview of supported and planned features.

#### Families of Standard Methods

- [x] Explicit Runge-Kutta Methods (ERK),
- [x] Explicit Partitioned Runge-Kutta Methods (EPRK),
- [x] Implicit Partitioned Runge-Kutta Methods (IPRK),
- [ ] Diagonally Implicit Runge-Kutta Methods (DIRK),
- [x] Fully Implicit Runge-Kutta Methods (FIRK),
- [ ] Singly Implicit Runge-Kutta Methods (SIRK),
- [ ] Additive Runge-Kutta Methods (ARK),
- [ ] Partitioned Additive Runge-Kutta Methods (PARK),
- [ ] Special Additive Runge-Kutta Methods (SARK),
- [ ] Special Partitioned Additive Runge-Kutta Methods (SPARK),
- [ ] Two-Step Runge-Kutta Methods (TSRK),
- [ ] General Linear Methods (GLM).

#### Families of Geometric Integrators

- [x] Gauss-Legendre Runge-Kutta Methods (GLRK),
- [x] Variational Partitioned Runge-Kutta Methods (VPRK),
- [ ] Variational Partitioned Additive Runge-Kutta Methods (VPARK),
- [ ] Continuous and Discontinuous Galerkin Variational Integrators (CGVI/DGVI),
- [ ] Spline Variational Integrators (SVI),
- [ ] Taylor Variational Integrators (TVI),
- [ ] Hamilton-Pontryagin-Galerkin Integrators (HPGI),
- [ ] Splitting Methods (SM).

#### Families of Equations

- [x] Systems of ODEs,
- [x] Systems of DAEs,
- [x] Partitioned ODEs,
- [x] Partitioned DAEs,
- [x] Implicit ODEs,
- [x] Implicit DAEs,

which can be prescribed manually or obtained as

- [ ] Euler-Lagrange Equations,
- [ ] Hamilton Equations,
- [ ] Hamilton-Pontryagin Equations,
- [ ] Lagrange-d'Alembert Equations,
- [ ] Hamilton-d'Alembert Equations,
- [ ] Symplectic Equations,
- [ ] Poisson Equations,

with

- [ ] Holonomic Constraints,
- [ ] Nonholonomic Constraints,
- [ ] Dirac Constraints.

#### Linear Solvers

- [x] LU decomposition (LAPACK),
- [x] LU decomposition (native Julia),
- [ ] Krylov,

#### Nonlinear Solvers

- [ ] Fixed-Point Iteration,
- [ ] Fixed-Point Iteration with Aitken's Acceleration,
- [ ] Jacobian-free Newton-Krylov,
- [x] Newton's method,
- [ ] Newton's method with line search (Armijo),
- [x] Quasi-Newton,

with

- [ ] Analytic Jacobian,
- [x] Finite Difference Jacobian,
- [ ] Jacobian obtained via Automatic Differentiation.

#### Diagnostics

- [ ] Runge-Kutta Stability Area,
- [ ] Convergence Analysis,
- [ ] First Poincaré Integral Invariant,
- [ ] Second Poincaré Integral Invariant.

#### Example Problems

- [x] Exponential Growth,
- [x] Pendulum,
- [x] Lotka-Volterra in 2D,
- [ ] Lotka-Volterra in 3D,
- [ ] Charged Particle in a Uniform Magnetic Field.


## License

The GeometricIntegrators.jl package is licensed under the [MIT "Expat" License](LICENSE.md).


# GeometricIntegrators.jl

*Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.*

[![Build Status](https://travis-ci.org/DDMGNI/GeometricIntegrators.jl.svg?branch=master)](https://travis-ci.org/DDMGNI/GeometricIntegrators.jl)
[![Coverage Status](https://coveralls.io/repos/github/DDMGNI/GeometricIntegrators.jl/badge.svg)](https://coveralls.io/github/DDMGNI/GeometricIntegrators.jl)
[![codecov](https://codecov.io/gh/DDMGNI/GeometricIntegrators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/DDMGNI/GeometricIntegrators.jl)

GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles.


## Manual

```@contents
Pages = ["tutorial.md",
         "integrators.md"]
```


## Modules

```@contents
Pages = ["modules/basis_functions.md",
         "modules/equations.md",
         "modules/integrators.md",
         "modules/interpolation.md",
         "modules/quadratures.md",
         "modules/discontinuities.md",
         "modules/simulations.md",
         "modules/solvers_linear.md",
         "modules/solvers_nonlinear.md",
         "modules/solutions.md",
         "modules/tableaus.md"
]
```


## Features

The following list provides and overview of supported and planned features.

#### Families of Standard Methods

- [x] Explicit Runge-Kutta Methods (ERK),
- [ ] Diagonally Implicit Runge-Kutta Methods (DIRK),
- [ ] Singly Diagonally Implicit Runge-Kutta Methods (SIRK),
- [x] Fully Implicit Runge-Kutta Methods (FIRK),
- [x] Explicit Partitioned Runge-Kutta Methods (EPRK),
- [x] Implicit Partitioned Runge-Kutta Methods (IPRK),
- [ ] Additive Runge-Kutta Methods (ARK),
- [ ] Specialised Additive Runge-Kutta Methods (SARK),
- [x] Partitioned Additive Runge-Kutta Methods (PARK),
- [ ] Specialised Partitioned Additive Runge-Kutta Methods (SPARK),
- [ ] Generalised Partitioned Additive Runge-Kutta Methods (GPARK),
- [ ] Two-step Runge-Kutta Methods (TSRK),
- [ ] General Linear Methods (GLM).

#### Families of Geometric Integrators

- [x] Gauss-Legendre Runge-Kutta Methods (GLRK),
- [x] Variational Partitioned Runge-Kutta Methods (VPRK),
- [x] Variational Partitioned Additive Runge-Kutta Methods (VPARK),
- [ ] Specialised Partitioned Additive Runge-Kutta Methods for Dirac Constraints (SPARK-Dirac),
- [x] Continuous Galerkin Variational Integrators (CGVI),
- [ ] Discontinuous Galerkin Variational Integrators (DGVI),
- [ ] Hamilton-Pontryagin-Galerkin Integrators (HPGI),
- [ ] Spline Variational Integrators (SVI),
- [ ] Taylor Variational Integrators (TVI),
- [x] Splitting Methods (SM).

#### Families of Stochastic Integrators

- [x] Stochastic Explicit Runge-Kutta Methods
- [x] Stochastic Implicit Runge-Kutta Methods
- [x] Stochastic Implicit Partitioned Runge-Kutta Methods
- [x] Stochastic Implicit Split Partitioned Runge-Kutta Methods
- [x] Stochastic Weak Explicit Runge-Kutta Methods
- [x] Stochastic Weak Implicit Runge-Kutta Methods

#### Families of Equations

- [x] Systems of ODEs,
- [x] Systems of DAEs,
- [x] Systems of SDEs,
- [x] Partitioned ODEs,
- [x] Partitioned DAEs,
- [ ] Partitioned SDEs,
- [x] Implicit ODEs,
- [x] Implicit DAEs,
- [ ] Implicit SDEs,
- [x] Hamiltonian DAEs,
- [x] Split ODEs,

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
- [ ] Fixed-Point Iteration with Aitken Acceleration,
- [ ] Fixed-Point Iteration with Anderson Acceleration,
- [ ] Jacobian-free Newton-Krylov,
- [x] Newton's method,
- [x] Newton's method with line search (Armijo, quadratic),
- [x] Quasi-Newton,

with

- [x] Finite Difference Jacobian,
- [x] Jacobian obtained via Automatic Differentiation.

#### Diagnostics

- [x] Symplecticity Conditions,
- [ ] Runge-Kutta Stability Area,
- [ ] Convergence Analysis,
- [x] First Poincaré Integral Invariant,
- [x] Second Poincaré Integral Invariant.


## References (mostly in preparation)

- Michael Kraus. Hamilton-Pontryagin-Galerkin Integrators.
- Michael Kraus. Projected Variational Integrators for Degenerate Lagrangian Systems.
- Michael Kraus. Variational Integrators for Noncanonical Hamiltonian Systems.
- Michael Kraus. Discontinuous Galerkin Variational Integrators for Degenerate Lagrangian Systems.
- Michael Kraus. Discontinuous Galerkin Variational Integrators for Hamiltonian Systems subject to Dirac Constraints.
- Michael Kraus. SPARK Methods for Degenerate Lagrangian Systems.
- Michael Kraus. SPARK Methods for Hamiltonian Systems subject to Dirac Constraints.
- Michael Kraus. Symplectic Runge-Kutta Methods for Certain Degenerate Lagrangian Systems.
- Michael Kraus and Tomasz M. Tyranowski􏰁. Variational Integrators for Stochastic Dissipative Hamiltonian Systems.


## Background Material

- Ernst Hairer and Christian Lubich. Numerical Solution of Ordinary Differential Equations. The Princeton Companion to Applied Mathematics, 293-305, 2015. Princeton University Press. ([Author's Web Site](https://na.uni-tuebingen.de/~lubich/pcam-ode.pdf))
- Ernst Hairer, Christian Lubich and Gerhard Wanner. Geometric Numerical Integration Illustrated by the Störmer–Verlet Method. Acta Numerica 12, 399-450, 2003. ([Journal](http://dx.doi.org/10.1017/S0962492902000144))
- Laurent O. Jay. Lobatto Methods. Encyclopedia of Applied and Computational Mathematics, 817–826. Springer, 2015. ([Article](http://dx.doi.org/10.1007/978-3-540-70529-1_123))


## Books on the Numerical Integration of Ordinary Differential Equations

- Ernst Hairer, Syvert P. Nørsett and Gerhard Wanner. Solving Ordinary Differential Equations I: Nonstiff Problems. Springer, 1993. ([eBook](http://link.springer.com/book/10.1007%2F978-3-540-78862-1))
- Ernst Hairer and Gerhard Wanner. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems. Springer, 1996. ([eBook](http://link.springer.com/book/10.1007%2F978-3-642-05221-7))
- Peter Deuflhard, Folkmar Bornemann. Scientific Computing with Ordinary Differential Equations. Springer, 2002. ([eBook](http://link.springer.com/book/10.1007/978-0-387-21582-2))
- John C. Butcher. Numerical Methods for Ordinary Differential Equations. Wiley, 2016. ([eBook](http://onlinelibrary.wiley.com/book/10.1002/9781119121534))
- Ernst Hairer, Christian Lubich and Gerhard Wanner. Geometric Numerical Integration. Springer, 2006. ([eBook](http://link.springer.com/book/10.1007%2F3-540-30666-8))
- Benedict Leimkuhler and Sebastian Reich. Simulating Hamiltonian Dynamics. Cambridge University Press, 2005. ([eBook](http://ebooks.cambridge.org/ebook.jsf?bid=CBO9780511614118))
- Sergio Blanes, Fernando Casas. A Concise Introduction to Geometric Numerical Integration. CRC Press, 2016. ([eBook](http://dx.doi.org/10.1201/b21563))


## License

> Copyright (c) 2016-2018 Michael Kraus <michael.kraus@ipp.mpg.de>
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

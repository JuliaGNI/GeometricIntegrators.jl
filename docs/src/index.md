
# GeometricIntegrators.jl

*Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.*

[![Build Status](https://travis-ci.org/DDMGNI/GeometricIntegrators.jl.svg?branch=master)](https://travis-ci.org/DDMGNI/GeometricIntegrators.jl)
[![Coverage Status](https://coveralls.io/repos/github/DDMGNI/GeometricIntegrators.jl/badge.svg)](https://coveralls.io/github/DDMGNI/GeometricIntegrators.jl)
![codecov](https://codecov.io/gh/DDMGNI/GeometricIntegrators.jl/branch/master/graph/badge.svg)

GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles.

## Features

Currently the following methods are supported or planned to be supported,

- [x] Explicit Runge-Kutta Methods (ERK),
- [x] Explicit Partitioned Runge-Kutta Methods (EPRK),
- [x] Implicit Partitioned Runge-Kutta Methods (IPRK),
- [x] Fully Implicit Runge-Kutta Methods (FIRK),
- [ ] Diagonally Implicit Runge-Kutta Methods (DIRK),
- [ ] Singly Implicit Runge-Kutta Methods (SIRK),
- [ ] Special Additive Runge-Kutta Methods (SARK),
- [ ] Special Partitioned Additive Runge-Kutta Methods (SPARK),
- [ ] Two-Step Runge-Kutta Methods (TSRK),
- [ ] General Linear Methods (GLM),
- [ ] Splitting Methods (SM).

The following families of equations are supported or planned to be supported,

- [x] Systems of ODEs,
- [ ] Systems of DAEs,
- [x] Partitioned ODEs,
- [ ] Partitioned DAEs,
- [x] Implicit ODEs
- [ ] Implicit DAEs

which can be prescribed manually or obtained as

- [ ] Euler-Lagrange Equations,
- [ ] Hamilton Equations,
- [ ] Symplectic Equations,
- [ ] Poisson Equations,

with

- [ ] Holonomic Constraints,
- [ ] Nonholonomic Constraints,
- [ ] Dirac Constraints.

The following families of integrators are supported or planned to be supported,

- [ ] Gauss-Legendre Runge-Kutta,
- [ ] Galerkin Variational Integrators,
- [ ] Taylor Variational Integrators,
- [ ] Hamilton-Pontryagin-Galerkin Integrators.

Available linear solvers are

- [x] LU decomposition (LAPACK),
- [ ] LU decomposition (native Julia),
- [ ] Krylov,

and nonlinear solvers

- [ ] Fixed-Point Iteration,
- [ ] Fixed-Point Iteration with Aitken's Acceleration,
- [x] Newton's method,
- [ ] Newton's method with line search,
- [x] Quasi-Newton,

either with exact Jacobian or with approximate Jacobian obtained via

- [x] Finite Differences,
- [ ] Automatic Differentiation,

and

- [ ] Jacobian-free Newton-Krylov.


## Manual

```@contents
Pages = ["tutorial.md"]
```


## Modules

```@contents
Pages = ["modules/equations.md",
         "modules/integrators.md",
         "modules/solvers_linear.md",
         "modules/solvers_nonlinear.md",
         "modules/tableaus.md"
]
```


## License

> Copyright (c) 2016 Michael Kraus <michael.kraus@ipp.mpg.de>
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


# GeomDAE.jl

GeomDAE.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeomDAE.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order e.g. via discrete variational principles.


## Features

Currently the following methods are supported or planned to be supported,

- [x] Explicit Runge-Kutta Methods (ERK),
- [x] Explicit Partitioned Runge-Kutta Methods (PRK),
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
- [ ] Discontinuous Galerkin Hamilton-Pontryagin Equations.

Available linear solvers are

- [x] LU decomposition (LAPACK),
- [ ] LU decomposition (native Julia),
- [ ] Krylov,

and nonlinear solvers

- [ ] Fixed-Point Iteration,
- [ ] Fixed-Point Iteration with Aitken's Acceleration,
- [x] Newton's method,
- [ ] Newton's method with line search,

either with exact Jacobian or with approximate Jacobian obtained via

- [x] Finite Differences,
- [ ] Automatic Differentiation,

and

- [ ] Jacobian-free Newton-Krylov.


## License

The GeomDAE.jl package is licensed under the MIT "Expat" License.

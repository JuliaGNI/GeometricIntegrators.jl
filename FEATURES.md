
# GeometricIntegrators.jl • Features

The following list provides and overview of supported and planned features.

#### Families of Standard Methods

- [x] Explicit Runge-Kutta Methods (ERK),
- [x] Diagonally Implicit Runge-Kutta Methods (DIRK),
- [x] Fully Implicit Runge-Kutta Methods (FIRK),
- [x] Gauss-Legendre Runge-Kutta Methods (GLRK),
- [x] Radau I, Radau II and Lobatto III Runge-Kutta Methods,
- [x] Explicit Partitioned Runge-Kutta Methods (EPRK),
- [x] Implicit Partitioned Runge-Kutta Methods (IPRK),
- [ ] Additive Runge-Kutta Methods (ARK),
- [x] Partitioned Additive Runge-Kutta Methods (PARK),
- [ ] Generalised Additive Runge-Kutta Methods (GARK),
- [x] Specialised Partitioned Additive Runge-Kutta Methods (SPARK),
- [ ] Continuous-stage Runge-Kutta Methods (CSRK),
- [ ] Two-step Runge-Kutta Methods (TSRK),
- [ ] General Linear Methods (GLM).

#### Families of Geometric Integrators

- [x] Variational Partitioned Runge-Kutta Methods (VPRK),
- [x] Hamiltonian Partitioned Additive Runge-Kutta Methods (HPARK, HSPARK),
- [x] Variational Partitioned Additive Runge-Kutta Methods (VPARK, VSPARK),
- [x] Continuous Galerkin Variational Integrators (CGVI),
- [x] Discontinuous Galerkin Variational Integrators (DGVI),
- [ ] Hamilton-Pontryagin-Galerkin Integrators (HPGI),
- [ ] Spline Variational Integrators (SVI),
- [ ] Taylor Variational Integrators (TVI),
- [x] Degenerate Variational Integrators (DVI),
- [x] Splitting Methods (SM),
- [ ] Hamiltonian Boundary Value Methods (HBVM).

#### Families of Stochastic Integrators

- [x] Stochastic Explicit Runge-Kutta Methods (SERK),
- [x] Stochastic Implicit Runge-Kutta Methods (SIRK),
- [x] Stochastic Implicit Partitioned Runge-Kutta Methods (SIPRK),
- [x] Stochastic Implicit Partitioned Additive Runge-Kutta Methods,
- [x] Stochastic Weak Explicit Runge-Kutta Methods (WERK),
- [x] Stochastic Weak Implicit Runge-Kutta Methods (WIRK).

#### Families of Equations

- [x] Systems of ODEs,
- [x] Systems of DAEs,
- [x] Systems of SDEs,
- [x] Partitioned ODEs,
- [x] Partitioned DAEs,
- [x] Partitioned SDEs,
- [x] Implicit ODEs,
- [x] Implicit DAEs,
- [ ] Implicit SDEs,
- [x] Variational ODEs,
- [x] Hamiltonian DAEs,
- [x] Split ODEs,
- [ ] Split Partitioned ODEs,
- [x] Split Partitioned SDEs.

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

- [x] Analytic Jacobian,
- [x] Finite Difference Jacobian,
- [x] Jacobian obtained via Automatic Differentiation.

#### Diagnostics

- [x] Symplecticity Conditions,
- [ ] Runge-Kutta Stability Area,
- [ ] Convergence Analysis,
- [x] First Poincaré Integral Invariant,
- [x] Second Poincaré Integral Invariant.


# Release Notes


## 0.10.0

### Breaking Changes

* Adapt Solution HDF5 interface to default Julia argument order and naming conventions
* Extract HDF5 functionality from Solutions into separate data structure
* Remove parallel Solution types

### New Features


### Documentation



## 0.9.0

### Breaking Changes

* Move `HermiteInterpolation` to Integrators and remove `Interpolation` sub-package
* Move `Equations` submodule to GeometricEquations.jl
* Move `Common`, `Config` and `Utils` submodules to GeometricBase.jl
* Move `TimeSeries`, `DataSeries` and `Solution` from `Solutions` types to GeometricBase.jl
* Remove parallel DataSeries and Solution types

### New Features

* Implement first and second order Degenerate Variational Integrators (DVIs)
* Add tests for extrapolation methods

### Fixes

* Bugfixes in implicit equations
* Bugfixes in extrapolation methods
* Bugfixes in initial guesses
* Bugfixes in VPRK and VSPARK initialisation
* Bugfixes in `TimeSeries` `getindex` methods

### Documentation

* Add missing docstrings in various places and remove superficial docstrings


## 0.8.0

### Breaking Changes

* Use RungeKutta.jl for most tableaus and coefficients
* Move stochastic integrators to separate package
* Rewrite of most equation types
* Rename `VODE` and `VDAE` to `LODE` and `LDAE` for consistency with `HODE` and `HDAE`
* Add optional fields for the secondary constraint to all *DAE equations

### New Features

* Allow for arbitrary data structures as states (still experimental and not fully supported)
* Add `convert` methods for `PODE` and `HODE` to `ODE` and `SODE`

### Fixes

* Countless minor bugfixes

### Documentation

* Add theoretical background for variational integrators, Runge-Kutta and splitting methods
* Add references for most methods


## 0.7.0

* Use CompactBasisFunctions.jl instead of BasisFunctions submodule
* Use QuadratureRules.jl instead of Quadratures submodule
* Use SimpleSolvers.jl instead of Solvers submodule
* Use GeometricProblems.jl instead of TestProblems submodule


## 0.6.2

* Bugfix release


## 0.6.1

* Bugfix release


## 0.6.0

### Breaking Changes

* Revise tableaus: align constructor names with RungeKutta.jl

### New Features

* Add new Runge-Kutta tableaus
* Generalise Lobatto and Radau tableaus to arbitrary number of stages
* Extend documentation on integrators and tableaus


## 0.5.1

* Update documentation
* Fix HDF5 v0.14 deprecations


## 0.5.0

* Moved repository to JuliaGNI
* Moved CI from Travis to GitHub

### Breaking Changes

* Functions for initial guesses are now called v̄ and f̄ and can be prescribed separately from v and f in PDAE, HDAE, etc.
* Rename SPARK tableau constructors and unify distinct constructors for Lobatto tableaus with different number of stages

### New Features

* Implement SPARK integrator for index-two DAEs
* Implement infrastructure for storing internal variables and solver output to atomic solutions
* Store internal variables of SPARK and VPRK integrators in atomic solution
* Add various five-stage Lobatto tableaus
* Add and clean up SPARK tableaus and add docstrings
* Add functions for checking symplecticity conditions of SPARK tableaus
* Add Aqua.jl tests

### Fixes

* Fix initial guess warnings in tests by prescribing proper functions for v̄ and f̄ in example problems
* Fix update_multiplier() method for SPARK integrators


## 0.4.1

### New Features

* Atomic solutions can now store a NamedTuple of internal variables of the integrator, including nonlinear solver output
* Output of internal variables has been added to VPRK integrators
* Add Gauss-Legendre tableaus for implicit partitioned Runge-Kutta methods

### Fixes

* Revision of integrator type hierarchy


## 0.4.0

### New Integrators

* Runge-Kutta integrators for implicit ODEs (`FIRKimplicit` and `SRKimplicit`)
* Variational Partitioned Runge-Kutta integrator with projection based on internal stages

### Fixes

* Computation of initial guess in *all* implicit integrators

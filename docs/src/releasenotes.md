
# Release Notes

## 0.4.1

### New Features

* Atomic solutions can now store a NamedTuple of internal variables of the integrator, including nonlinear solver output
* Output of internal variables has been added to VPRK integrators

### Fixes

* Revision of integrator type hierarchy


## 0.4.0

### New Integrators

* Runge-Kutta integrators for implicit ODEs (`FIRKimplicit` and `SRKimplicit`)
* Variational Partitioned Runge-Kutta integrator with projection based on internal stages

### Fixes

* Computation of initial guess in *all* implicit integrators

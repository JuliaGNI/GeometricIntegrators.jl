# Adaptive Time Stepping

GeometricIntegrators.jl does not provide any general infrastructure for adaptive time stepping.
The main reason is that adaptive time stepping is not easy to combine with structure-preservation.
Most typical applicatons of the GeometricIntegrators developers require the output of solutions at a sequence of time steps with constant step size.

Nonetheless, geometric methods with time step adaptation do exist and it is not hard to implement such methods in the GeometricIntegrators framework.
Here, the standard infrastructure is used to specify "target time steps", at which a solution has to be computed.
That means an adaptive integrator is required to compute a solution for every point in the equidistant time series, but in between in adapts as it wishes. Whenever an adaptive time step would step over a "target time step", it is reduced to hit that target.
For geometric, structure-preserving integrators that often is the only sensible thing to do (and for most practical applications as well).

In all of this, the integrator has to take care of the sub-cycling and the rest of GeometricIntegrators doesn't really care about it. If for some reason you want to output the solution at the intermediate irregular time steps, this is relatively easily possible via the (still mostly undocumented) [mid-level interface](code_integration.md) that is used to call the integrators.

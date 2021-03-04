# Equations

In *GeometricIntegrators.jl* we support three basic types of equations:
* ordinary differential equations (ODEs),
* differential algebraic equations (DAEs),
* stochastic differential equations (SDEs).

For each type, there are several subtypes
* standard equations ([`ODE`](@ref), [`DAE`](@ref), [`SDE`](@ref)),
* implicit equations ([`IODE`](@ref), [`IDAE`](@ref)),
* partitioned equations ([`PODE`](@ref), [`PDAE`](@ref), [`PSDE`](@ref)),
* Hamiltonian equations ([`HODE`](@ref), [`HDAE`](@ref)),
* Lagrangian equations ([`LODE`](@ref), [`LDAE`](@ref)),
* split equations ([`SODE`](@ref), [`SPDAE`](@ref)), [`SPSDE`](@ref)).

Each equation holds a number of functions determining the vector field, constraints, initial conditions, and possibly additional information like parameters, periodicity, invariants and the Hamiltonian or Lagrangian.

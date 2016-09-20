
# GeomDAE.jl

GeomDAE.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeomDAE.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order e.g. via discrete variational principles.



## License

The GeomDAE.jl package is licensed under the MIT "Expat" License.

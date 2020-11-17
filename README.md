
# GeometricIntegrators.jl

*Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.*

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/GeometricIntegrators.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliagni.github.io/GeometricIntegrators.jl/latest/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
![CI](https://github.com/JuliaGNI/GeometricIntegrators.jl/workflows/CI/badge.svg)
[![Build Status](https://travis-ci.org/JuliaGNI/GeometricIntegrators.jl.svg?branch=master)](https://travis-ci.org/JuliaGNI/GeometricIntegrators.jl)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricIntegrators.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricIntegrators.html)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNI/GeometricIntegrators.jl/badge.svg)](https://coveralls.io/github/JuliaGNI/GeometricIntegrators.jl)
[![codecov Status](https://codecov.io/gh/JuliaGNI/GeometricIntegrators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/GeometricIntegrators.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3648325.svg)](https://doi.org/10.5281/zenodo.3648325)

GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles.  

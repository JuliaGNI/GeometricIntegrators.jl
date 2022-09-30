---
title: 'GeometricIntegrators.jl: Geometric Numerical Integration in Julia'
tags:
  - Julia
  - differential equations
  - Lagrangian dynamics
  - Hamiltonian dynamics
  - symplectic integrators
  - variational integrators
authors:
  - name: Michael Kraus
    orcid: 0000-0001-7637-3979
    affiliation: "1, 2"
affiliations:
 - name: Max-Planck-Institut für Plasmaphysik, Boltzmannstraße 2, 85748 Garching, Germany
   index: 1
 - name: Technische Universität München, Zentrum Mathematik, Boltzmannstraße 3, 85748 Garching, Germany
   index: 2
date: 23 November 2020
bibliography: paper.bib
---

# Summary

*GeometricIntegrators.jl* is a library of geometric numerical algorithms for the solution of differential equations in Julia.
It provides a unified interface for many different algorithms and various types of equations, such as ordinary differential equations, stochastic differential equations and differential algebraic equations, in particular in Lagrangian (variational) and Hamiltonian (symplectic) form.
It aims at providing a comprehensive collection of geometric or structure-preserving algorithms, which can be used either interactively or as computational core in other codes.
The library provides both, a high-level interface that requires only very few lines of code to solve an actual problem, and a lean low-level interface that allows for straightforward integration into application codes via the exchange of very small data structures.
In both cases, the library leaves maximum control to the user, e.g., with respect to the choice of numerical methods and the setup of linear and nonlinear solvers.


# Statement of need

Differential equations are ubiquitous in science and engineering. Many equations possess geometric features or abstract mathematical structures that need to be preserved in the discretisation in order to obtain reliable simulation results, especially for nonlinear problems and long-time simulations. The preservation of such properties improves stability, bounds global error growth and reduces numerical artefacts [@BlanesCasas:2016; @HairerLubichWanner:2006; @LeimkuhlerReich:2004; @SanzSernaCalvo:1994].
Robust, performant and structure-preserving solvers for different types of differential equations are thus needed across many disciplines. *GeometricIntegrators.jl* provides such solvers and makes them available for both direct use as well as integration into other codes. Furthermore, the implemented algorithms can also be used within the *DifferentialEquations.jl* ecosystem [@Rackauckas:2017], which is the defacto standard differential equation solver for the Julia programming language [@Bezanson:2017].

*GeometricIntegrators.jl* provides a comprehensive library of geometric integration algorithms as well as some non-geometric algorithms.
It collects native Julia implementations of many known methods under a unified interface.
Once a problem is implemented in the GeometricIntegators framework, all of its algorithms can immediately be applied and their performance evaluated.
This facilitates numerical experiments with a wide variety of algorithms, simplifies benchmarking, and makes it easy to identify the best algorithm for a given problem.
The implemented algorithms include explicit, implicit and partitioned Runge-Kutta methods, SPARK methods, splitting methods, symplectic methods and variational integrators.
Most methods are implemented in an abstract way that allows for the flexible choice of tableaus, approximation spaces, basis functions, quadrature rules, and thus order of convergence. It also means that adding e.g. a new Runge-Kutta or splitting method merely amounts to adding its tableau.
As most geometric integrators are not easily combined with time step adaption in a structure-preserving way, GeometricIntegrators.jl does not provide any general infrastructure for adaptive time stepping. Nonetheless, individual integrators can implement their own adaptivity strategies as long as they provide a solution at a predefined, equidistant series of time steps.

One of the original aims of *GeometricIntegrators.jl* is to serve as a testbed for the development and analysis of novel algorithms [@Kraus:2017; @Kraus:2020; @KrausTyranowski:2020].
Due to the modular structure and the use of the multiple dispatch paradigm, the library can easily be extended, e.g., towards new algorithms or new types of equations.
The library aims at providing efficient implementations of diverse algorithms in order to be able to perform simulations and benchmarks with millions or even billions of time steps that facilitate the study of the long-time behaviour of both numerical algorithms and dynamical systems.
The current scope of applications is mainly small- to mid-size systems of differential equations, e.g., systems of ordinary differential equations or semidiscretisations of partial differential equations with a moderate number of degrees of freedom.
It is envisaged that in the future GeometricIntegrators.jl will also be able to address larger problems, especially semidiscretisations of partial differential equations in higher dimensions.
In particular, this requires interfaces to appropriate iterative and parallel linear solvers, which are still lacking.


# Other Software

A Julia package closely related to *GeometricIntegrators.jl* is *DifferentialEquations.jl* [@Rackauckas:2017]. Althouh the two libraries provide similar functionality, their scope is rather different. While *DifferentialEquations.jl* provides a feature-rich ecosystem for the solution of differential equations, the focus of *GeometricIntegrators.jl* is on algorithms. 
While *DifferentialEquations.jl* is based around adaptive time stepping algorithms, GeometricIntegrators.jl focuses on fixed time step methods, given that structure-preservation is not easily combined with time step adaptivity. Therefore both packages are rather complementary, and *GeometricIntegrators.jl* can in fact be used as backend for *DifferentialEquations.jl* via [GeometricIntegratorsDiffEq.jl](https://github.com/JuliaDiffEq/GeometricIntegratorsDiffEq.jl).


# Quality control and contributions

GeometricIntegrators.jl uses continuous integration testing with all supported versions of Julia on macOS, Linux and Windows via GitHub CI.
The tests cover most of the library, checking for the correct functioning of all submodules and convergence of the implemented algorithms.
The tests can also be run locally in the Julia's REPL by the command `]test GeometricIntegrators`.

Support and submission of contributions to the library are handled through the GitHub repository via issues or by pull requests.


# Acknowledgements

This work has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 708124. The views and opinions expressed herein do not necessarily reflect those of the European Commission.

Tomasz M. Tyranowski provided the initial implementation of stochastic integrators which can now be found in the package StochasticIntegrators.jl.

# References

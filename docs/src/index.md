
# GeometricIntegrators.jl

*Julia library of geometric integrators for differential equations.*

![CI](https://github.com/JuliaGNI/GeometricIntegrators.jl/workflows/CI/badge.svg)
[![Build Status](https://travis-ci.org/JuliaGNI/GeometricIntegrators.jl.svg?branch=master)](https://travis-ci.org/JuliaGNI/GeometricIntegrators.jl)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricIntegrators.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricIntegrators.html)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNI/GeometricIntegrators.jl/badge.svg)](https://coveralls.io/github/JuliaGNI/GeometricIntegrators.jl)
[![codecov Status](https://codecov.io/gh/JuliaGNI/GeometricIntegrators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/GeometricIntegrators.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3648325.svg)](https://doi.org/10.5281/zenodo.3648325)


GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations, stochastic differential equations and differential algebraic equations in Julia.
Its main purpose is the democratization and proliferation of geometric integrators by providing a comprehensive collection of structure-preserving as well as standard algorithms under a unified interface. 
GeometricIntegrators.jl can be used either interactively or as computational core in other codes. It provides both, a high-level interface that requires only very few lines of code to solve an actual problem, and a lean low-level interface that allows for straightforward integration into application codes via the exchange of minimalistic data structures.
In both, the library leaves maximum control to the user. While trying to pick sensible defaults, all settings are accessible to and modifiable by the user. Suitable abstraction layers allow to choose between different linear and nonlinear solvers, auto-differentiation packages or custom routines for the computation of Jacobians and the like.


## Statement of need

Differential equations are ubiquitous in science and engineering. Many equations possess geometric features or abstract mathematical structures that need to be preserved in the discretisation in order to obtain reliable simulation results, especially for nonlinear problems and long-time simulations. The preservation of such properties improves stability, bounds global error growth and reduces numerical artefacts.
Robust, performant and structure-preserving solvers for different types of differential equations are thus needed across many disciplines. *GeometricIntegrators.jl* provides such solvers and makes them available for both direct use as well as integration into other codes. Furthermore, the implemented algorithms can also be used within the *[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)* ecosystem [[Rackauckas:2017](@cite))], which is the defacto standard differential equation solver for the Julia programming language [[Bezanson:2017](@cite)].

*GeometricIntegrators.jl* provides a comprehensive library of existing geometric integration as well as non-geometric algorithms, such as explicit, implicit, partitioned and stochastic Runge-Kutta methods, SPARK methods, splitting methods, symplectic methods and variational integrators. Most methods are implemented in an abstract way that allows for the flexible choice of tableaus, approximation spaces, basis functions, quadrature rules, and thus order of convergence.
*GeometricIntegrators.jl* also serves as a testbed for the development and analysis of novel algorithms. Due to the modular structure and the use of the multiple dispatch paradigm, the library can easily be extended, e.g., towards new algorithms or new types of equations. The library is designed to minimize overhead and maximize performance in order to be able to perform simulations with millions or even billions of time steps to facilitate the study of the long-time behaviour of both numerical algorithms and dynamical systems.


## Other Software

A Julia package closely related to *GeometricIntegrators.jl* is *[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)* [[Rackauckas:2017](@cite)]. However, the scope of the two libraries is rather different. While *DifferentialEquations.jl* provides a feature-rich ecosystem for the solution of differential equations, the focus of *GeometricIntegrators.jl* is on algorithms. In fact, *GeometricIntegrators.jl* can be used as backend for *DifferentialEquations.jl* via [GeometricIntegratorsDiffEq.jl](https://github.com/JuliaDiffEq/GeometricIntegratorsDiffEq.jl).


## Manual

```@contents
Pages = ["tutorial/tutorial.md",
         "integrators.md"]
```


## Modules

```@contents
Pages = ["modules/basis_functions.md",
         "modules/equations.md",
         "modules/integrators.md",
         "modules/interpolation.md",
         "modules/quadratures.md",
         "modules/discontinuities.md",
         "modules/simulations.md",
         "modules/solvers_linear.md",
         "modules/solvers_nonlinear.md",
         "modules/solutions.md",
         "modules/tableaus.md"
]
```


## References

- Michael Kraus. Hamilton-Pontryagin-Galerkin Integrators.
- Michael Kraus. Projected Variational Integrators for Degenerate Lagrangian Systems. [arXiv:1708.07356](https://arxiv.org/abs/1708.07356).
- Michael Kraus. SPARK Methods for Degenerate Lagrangian Systems.
- Michael Kraus. SPARK Methods for Hamiltonian Systems subject to Dirac Constraints.
- Michael Kraus. Symplectic Runge-Kutta Methods for Degenerate Lagrangian Systems.
- Michael Kraus. Variational Integrators for Noncanonical Hamiltonian Systems.
- Michael Kraus and Tomasz M. Tyranowski. Variational Integrators for Stochastic Dissipative Hamiltonian Systems. [arXiv:1909.07202](https://arxiv.org/abs/1909.07202),
  [Journal](https://doi.org/10.1088/1742-6596/1391/1/012037).
- Michael Kraus, Joshua Burby. Conservation of Poincaré Integral Invariants in Numerical Simulations.


## Background Material

- Ernst Hairer and Christian Lubich. Numerical Solution of Ordinary Differential Equations. The Princeton Companion to Applied Mathematics, 293-305, 2015. Princeton University Press. ([Author's Web Site](https://na.uni-tuebingen.de/~lubich/pcam-ode.pdf))
- Ernst Hairer, Christian Lubich and Gerhard Wanner. Geometric Numerical Integration Illustrated by the Störmer–Verlet Method. Acta Numerica 12, 399-450, 2003. ([Journal](http://dx.doi.org/10.1017/S0962492902000144))
- Laurent O. Jay. Lobatto Methods. Encyclopedia of Applied and Computational Mathematics, 817–826. Springer, 2015. ([Article](http://dx.doi.org/10.1007/978-3-540-70529-1_123))


## Books on Geometric Numerical Integration

- Sergio Blanes, Fernando Casas. A Concise Introduction to Geometric Numerical Integration. CRC Press, 2016. ([eBook](http://dx.doi.org/10.1201/b21563))
- Ernst Hairer, Christian Lubich and Gerhard Wanner. Geometric Numerical Integration. Springer, 2006. ([eBook](http://link.springer.com/book/10.1007%2F3-540-30666-8))
- Benedict Leimkuhler and Sebastian Reich. Simulating Hamiltonian Dynamics. Cambridge University Press, 2005. ([eBook](http://ebooks.cambridge.org/ebook.jsf?bid=CBO9780511614118))
- Jesús Maria Sanz-Serna, Manuel P. Calvo. Numerical Hamiltonian Problems. Chapman Hall, 1994.


## Books on the Numerical Integration of Differential Equations

- Ernst Hairer, Syvert P. Nørsett and Gerhard Wanner. Solving Ordinary Differential Equations I: Nonstiff Problems. Springer, 1993. ([eBook](http://link.springer.com/book/10.1007%2F978-3-540-78862-1))
- Ernst Hairer and Gerhard Wanner. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems. Springer, 1996. ([eBook](http://link.springer.com/book/10.1007%2F978-3-642-05221-7))
- Ernst Hairer, Christian Lubich, Michel Roche. The Numerical Solution of Differential-Algebraic Systems by Runge-Kutta Methods. Springer, 1989. ([eBook](https://link.springer.com/book/10.1007/BFb0093947))
- Peter Deuflhard, Folkmar Bornemann. Scientific Computing with Ordinary Differential Equations. Springer, 2002. ([eBook](http://link.springer.com/book/10.1007/978-0-387-21582-2))
- John C. Butcher. Numerical Methods for Ordinary Differential Equations. Wiley, 2016. ([eBook](http://onlinelibrary.wiley.com/book/10.1002/9781119121534))


## License

> Copyright (c) Michael Kraus <michael.kraus@ipp.mpg.de>
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

var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#GeometricIntegrators.jl-1",
    "page": "Home",
    "title": "GeometricIntegrators.jl",
    "category": "section",
    "text": "Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.(Image: Build Status) (Image: Coverage Status) (Image: codecov)GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles."
},

{
    "location": "#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\"tutorial.md\",\n         \"integrators.md\"]"
},

{
    "location": "#Modules-1",
    "page": "Home",
    "title": "Modules",
    "category": "section",
    "text": "Pages = [\"modules/basis_functions.md\",\n         \"modules/equations.md\",\n         \"modules/integrators.md\",\n         \"modules/interpolation.md\",\n         \"modules/quadratures.md\",\n         \"modules/discontinuities.md\",\n         \"modules/simulations.md\",\n         \"modules/solvers_linear.md\",\n         \"modules/solvers_nonlinear.md\",\n         \"modules/solutions.md\",\n         \"modules/tableaus.md\"\n]"
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "The following list provides and overview of supported and planned features."
},

{
    "location": "#Families-of-Standard-Methods-1",
    "page": "Home",
    "title": "Families of Standard Methods",
    "category": "section",
    "text": "[x] Explicit Runge-Kutta Methods (ERK),\n[x] Diagonally Implicit Runge-Kutta Methods (DIRK),\n[x] Fully Implicit Runge-Kutta Methods (FIRK),\n[x] Explicit Partitioned Runge-Kutta Methods (EPRK),\n[x] Implicit Partitioned Runge-Kutta Methods (IPRK),\n[ ] Additive Runge-Kutta Methods (ARK),\n[ ] Specialised Additive Runge-Kutta Methods (SARK),\n[x] Partitioned Additive Runge-Kutta Methods (PARK),\n[ ] Specialised Partitioned Additive Runge-Kutta Methods (SPARK),\n[ ] Generalised Partitioned Additive Runge-Kutta Methods (GPARK),\n[ ] Two-step Runge-Kutta Methods (TSRK),\n[ ] General Linear Methods (GLM)."
},

{
    "location": "#Families-of-Geometric-Integrators-1",
    "page": "Home",
    "title": "Families of Geometric Integrators",
    "category": "section",
    "text": "[x] Gauss-Legendre Runge-Kutta Methods (GLRK),\n[x] Variational Partitioned Runge-Kutta Methods (VPRK),\n[x] Variational Partitioned Additive Runge-Kutta Methods (VPARK),\n[ ] Hamiltonian Partitioned Additive Runge-Kutta Methods (HPARK),\n[x] Continuous Galerkin Variational Integrators (CGVI),\n[x] Discontinuous Galerkin Variational Integrators (DGVI),\n[ ] Hamilton-Pontryagin-Galerkin Integrators (HPGI),\n[ ] Spline Variational Integrators (SVI),\n[ ] Taylor Variational Integrators (TVI),\n[x] Splitting Methods (SM),\n[ ] Hamiltonian Boundary Value Methods (HBVM)."
},

{
    "location": "#Families-of-Stochastic-Integrators-1",
    "page": "Home",
    "title": "Families of Stochastic Integrators",
    "category": "section",
    "text": "[x] Stochastic Explicit Runge-Kutta Methods,\n[x] Stochastic Implicit Runge-Kutta Methods,\n[x] Stochastic Implicit Partitioned Runge-Kutta Methods,\n[x] Stochastic Implicit Split Partitioned Runge-Kutta Methods,\n[x] Stochastic Weak Explicit Runge-Kutta Methods,\n[x] Stochastic Weak Implicit Runge-Kutta Methods."
},

{
    "location": "#Families-of-Equations-1",
    "page": "Home",
    "title": "Families of Equations",
    "category": "section",
    "text": "[x] Systems of ODEs,\n[x] Systems of DAEs,\n[x] Systems of SDEs,\n[x] Partitioned ODEs,\n[x] Partitioned DAEs,\n[x] Partitioned SDEs,\n[x] Implicit ODEs,\n[x] Implicit DAEs,\n[ ] Implicit SDEs,\n[x] Variational ODEs,\n[x] Hamiltonian DAEs,\n[x] Split ODEs,\n[x] Split Partitioned SDEs,which can be prescribed manually or obtained as[ ] Euler-Lagrange Equations,\n[ ] Hamilton Equations,\n[ ] Hamilton-Pontryagin Equations,\n[ ] Lagrange-d\'Alembert Equations,\n[ ] Hamilton-d\'Alembert Equations,\n[ ] Symplectic Equations,\n[ ] Poisson Equations,with[ ] Holonomic Constraints,\n[ ] Nonholonomic Constraints,\n[ ] Dirac Constraints."
},

{
    "location": "#Linear-Solvers-1",
    "page": "Home",
    "title": "Linear Solvers",
    "category": "section",
    "text": "[x] LU decomposition (LAPACK),\n[x] LU decomposition (native Julia),\n[ ] Krylov,"
},

{
    "location": "#Nonlinear-Solvers-1",
    "page": "Home",
    "title": "Nonlinear Solvers",
    "category": "section",
    "text": "[ ] Fixed-Point Iteration,\n[ ] Fixed-Point Iteration with Aitken Acceleration,\n[ ] Fixed-Point Iteration with Anderson Acceleration,\n[ ] Jacobian-free Newton-Krylov,\n[x] Newton\'s method,\n[x] Newton\'s method with line search (Armijo, quadratic),\n[x] Quasi-Newton,with[x] Analytic Jacobian,\n[x] Finite Difference Jacobian,\n[x] Jacobian obtained via Automatic Differentiation."
},

{
    "location": "#Diagnostics-1",
    "page": "Home",
    "title": "Diagnostics",
    "category": "section",
    "text": "[x] Symplecticity Conditions,\n[ ] Runge-Kutta Stability Area,\n[ ] Convergence Analysis,\n[x] First Poincaré Integral Invariant,\n[x] Second Poincaré Integral Invariant."
},

{
    "location": "#References-(mostly-in-preparation)-1",
    "page": "Home",
    "title": "References (mostly in preparation)",
    "category": "section",
    "text": "Michael Kraus. Hamilton-Pontryagin-Galerkin Integrators.\nMichael Kraus. Projected Variational Integrators for Degenerate Lagrangian Systems.\nMichael Kraus. Variational Integrators for Noncanonical Hamiltonian Systems.\nMichael Kraus. Discontinuous Galerkin Variational Integrators for Degenerate Lagrangian Systems.\nMichael Kraus. Discontinuous Galerkin Variational Integrators for Hamiltonian Systems subject to Dirac Constraints.\nMichael Kraus. SPARK Methods for Degenerate Lagrangian Systems.\nMichael Kraus. SPARK Methods for Hamiltonian Systems subject to Dirac Constraints.\nMichael Kraus. Symplectic Runge-Kutta Methods for Certain Degenerate Lagrangian Systems.\nMichael Kraus and Tomasz M. Tyranowski􏰁. Variational Integrators for Stochastic Dissipative Hamiltonian Systems."
},

{
    "location": "#Background-Material-1",
    "page": "Home",
    "title": "Background Material",
    "category": "section",
    "text": "Ernst Hairer and Christian Lubich. Numerical Solution of Ordinary Differential Equations. The Princeton Companion to Applied Mathematics, 293-305, 2015. Princeton University Press. (Author\'s Web Site)\nErnst Hairer, Christian Lubich and Gerhard Wanner. Geometric Numerical Integration Illustrated by the Störmer–Verlet Method. Acta Numerica 12, 399-450, 2003. (Journal)\nLaurent O. Jay. Lobatto Methods. Encyclopedia of Applied and Computational Mathematics, 817–826. Springer, 2015. (Article)"
},

{
    "location": "#Books-on-the-Numerical-Integration-of-Ordinary-Differential-Equations-1",
    "page": "Home",
    "title": "Books on the Numerical Integration of Ordinary Differential Equations",
    "category": "section",
    "text": "Ernst Hairer, Syvert P. Nørsett and Gerhard Wanner. Solving Ordinary Differential Equations I: Nonstiff Problems. Springer, 1993. (eBook)\nErnst Hairer and Gerhard Wanner. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems. Springer, 1996. (eBook)\nPeter Deuflhard, Folkmar Bornemann. Scientific Computing with Ordinary Differential Equations. Springer, 2002. (eBook)\nJohn C. Butcher. Numerical Methods for Ordinary Differential Equations. Wiley, 2016. (eBook)\nErnst Hairer, Christian Lubich and Gerhard Wanner. Geometric Numerical Integration. Springer, 2006. (eBook)\nBenedict Leimkuhler and Sebastian Reich. Simulating Hamiltonian Dynamics. Cambridge University Press, 2005. (eBook)\nSergio Blanes, Fernando Casas. A Concise Introduction to Geometric Numerical Integration. CRC Press, 2016. (eBook)"
},

{
    "location": "#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "Copyright (c) 2016-2018 Michael Kraus <michael.kraus@ipp.mpg.de>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

{
    "location": "tutorial/#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": "In the simplest cases, the use of GeometricIntegrators.jl requires the construction of two objects, an equation and an integrator. The integrator is usually implicitly selected by specifying an equation and a tableau."
},

{
    "location": "tutorial/#Equations-1",
    "page": "Tutorial",
    "title": "Equations",
    "category": "section",
    "text": "In GeometricIntegrators.jl we distinguish between three basic types of equations: ordinary differential equations (ODEs), differential algebraic equations (DAEs) and stochastic differential equations (SDEs). For each type, there are several subtypes like implicit equations (IODE, etc.), partitioned equations (PODE, etc.) or split equations (SODE, etc.).Instantiating an ODE object for the pendulum problem \\[ \\dot{x}1 = x2 , \\hspace{3em} \\dot{x}2 = \\sin (x1) , \\] can be achieved byfunction pendulum_rhs(t, x, f)\n    f[1] = x[2]\n    f[2] = sin(x[1])\nend\n\node = ODE(pendulum_rhs, [acos(0.4), 0.0])The first argument to the ODE constructor is the function that determines the vector field of the equation dotx (t) = f(t x(t)), and the second argument determines the initial conditions. The function defining the vector field has to take three arguments, the current time t, the current solution vector x and the output vector f.The pendulum problem is a Hamiltonian system that can also be expressed as \\[ \\dot{q} = \\frac{\\partial H}{\\partial p} = p , \\hspace{3em} \\dot{p} = - \\frac{\\partial H}{\\partial q} = \\sin (q) , \\hspace{3em} H (q,p) = \\frac{1}{2} p^2 + \\cos (q) . \\] This structure, namely the partitioning into two sets of variables (qp) instead of x, can be exploited for more efficient integration. Such equations can be defined in terms of a partitioned ODE, where the vector fields are specified separately,function pendulum_v(t, q, p, v)\n    v[1] = p[1]\nend\n\nfunction pendulum_f(t, q, p, f)\n    f[1] = sin(q[1])\nend\n\npode = PODE(pendulum_v, pendulum_f, [acos(0.4)], [0.0])The first two arguments to the PODE constructor are the functions that determine the vector fields of the equations dotq (t) = v(t q(t) p(t)) and dotp (t) = f(t q(t) p(t)). The third and fourth argument determines the initial conditions of q and p, respectively. The functions defining the vector field have to take four arguments, the current time t, the current solution vectors q and p and the output vector v or f."
},

{
    "location": "tutorial/#Integrators-1",
    "page": "Tutorial",
    "title": "Integrators",
    "category": "section",
    "text": "We support a number of standard integrators (geometric and non-geometric) like explicit, implicit and partitioned Runge-Kutta methods, splitting methods and general linear methods (planned).In order to instantiate many of the standard integrators, one needs to specify an ODE, a tableau and a timestep, e.g.,int = Integrator(ode, getTableauExplicitEuler(), 0.1)In order to run the integrator, the integrate() functions is called, passing an integrator object and the number of time steps to integrate:sol = integrate(int, 10)The integrate function automatically creates an appropriate solution object, that contains the result of the integration.For a Hamiltonian system, defined as a PODE, a different tableau might be more appropriate, for example a symplectic Euler method,int = Integrator(pode, getTableauSymplecticEulerA(), 0.1)\nsol = integrate(int, 10)This creates a different integrator, which exploits the partitioned structure of the system. The solution return by the integrate step will also be a different solution, adapted to the partitioned system."
},

{
    "location": "tutorial/#Tableaus-1",
    "page": "Tutorial",
    "title": "Tableaus",
    "category": "section",
    "text": "Many tableaus for Runge-Kutta methods are predefined and can easily be used like outlined above. In particular, this includes the following methods:"
},

{
    "location": "tutorial/#Explicit-Runge-Kutta-Methods-1",
    "page": "Tutorial",
    "title": "Explicit Runge-Kutta Methods",
    "category": "section",
    "text": "Function Order Method\ngetTableauExplicitEuler() 1 Explicit / Forward Euler\ngetTableauExplicitMidpoint() 2 Explicit Midpoint\ngetTableauHeun() 2 Heun\'s Method\ngetTableauKutta() 3 Kutta\'s Method\ngetTableauERK4() 4 Explicit 4th order Runge-Kutta (1/6 rule)\ngetTableauERK438() 4 Explicit 4th order Runge-Kutta (3/8 rule)"
},

{
    "location": "tutorial/#Fully-Implicit-Runge-Kutta-Methods-1",
    "page": "Tutorial",
    "title": "Fully Implicit Runge-Kutta Methods",
    "category": "section",
    "text": "Function Order Method\ngetTableauImplicitEuler() 1 Implicit / Backward Euler\ngetTableauImplicitMidpoint() 2 Implicit Midpoint\ngetTableauRadIIA2() 3 Radau-IIA s=2\ngetTableauRadIIA3() 5 Radau-IIA s=3\ngetTableauSRK3() 4 Symmetric Runge-Kutta s=3\ngetTableauGLRK(s) 2s Gauss-Legendre Runge-Kutta"
},

{
    "location": "tutorial/#Explicit-Partitioned-Runge-Kutta-Methods-1",
    "page": "Tutorial",
    "title": "Explicit Partitioned Runge-Kutta Methods",
    "category": "section",
    "text": "Function Order Method\ngetTableauSymplecticEulerA() 1 Symplectic Euler A\ngetTableauSymplecticEulerB() 1 Symplectic Euler B\ngetTableauLobattoIIIAIIIB2() 2 Lobatto-IIIA-IIIB\ngetTableauLobattoIIIBIIIA2() 2 Lobatto-IIIB-IIIA"
},

{
    "location": "tutorial/#Custom-Tableaus-1",
    "page": "Tutorial",
    "title": "Custom Tableaus",
    "category": "section",
    "text": "If required, it is straight-forward to create a custom tableau. The tableau of Heun\'s method, for example, is defined as follows:a = [[0.0 0.0]\n     [1.0 0.0]]\nb = [0.5, 0.5]\nc = [0.0, 1.0]\no = 2\n\ntab = TableauERK(:heun, o, a, b, c)Here, o is the order of the method, a are the coefficients, b the weights and c the nodes. TableauERK states that the method is explicit. Other choices include TableauFIRK for fully implicit Runge-Kutta methods, TableauDIRK for diagonally implicit and TableauSIRK for singly implicit Runge-Kutta methods. TableauEPRK and TableauIPRK can be used for explicit and implicit partitioned Runge-Kutta methods. The first parameter of the constructor of each tableau assigns a name to the tableau. Such custom tableaus can be used in exactly the same as standard tableaus, e.g., byint = Integrator(ode, tab, 0.1)\nsol = integrate(int, 10)making it very easy to implement and test new methods."
},

{
    "location": "tutorial/#Solutions-1",
    "page": "Tutorial",
    "title": "Solutions",
    "category": "section",
    "text": "In what we have seen so far, the solution was always automatically created by the integrate() function. While this is often convenient, it is sometimes not performant, e.g., when carrying out long-time simulations with intermediate saving of the solution. In such cases, it is better to preallocate a solution object bysol = Solution(ode, 0.1, 10)where the first argument is an equation, the second argument is the time step and the third argument is the number of time steps that will be computed in one integration step. The call to the integrator is then made viaintegrate!(int, sol)If several integration cycles shall be performed, the reset!() function can be used to copy the solution of the last time step to the initial conditions of the solution,for i in 1:10\n    integrate!(int, sol)\n    #\n    # save or process solution\n    #\n    reset!(sol)\nendAll solutions have a t field holding the series of time steps that has been computed in addition to several data fields, for example q for an ODE solution, q and p for a PODE solution, qand λ for a DAE solution, and q, p and λ for a PDAE solution."
},

{
    "location": "integrators/#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/#Integrators-1",
    "page": "Overview",
    "title": "Integrators",
    "category": "section",
    "text": ""
},

{
    "location": "integrators/splitting/#",
    "page": "Splitting",
    "title": "Splitting",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/splitting/#Splitting-Methods-1",
    "page": "Splitting",
    "title": "Splitting Methods",
    "category": "section",
    "text": ""
},

{
    "location": "integrators/rk/#",
    "page": "Runge-Kutta",
    "title": "Runge-Kutta",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/rk/#Runge-Kutta-Methods-1",
    "page": "Runge-Kutta",
    "title": "Runge-Kutta Methods",
    "category": "section",
    "text": ""
},

{
    "location": "integrators/rk/#Gauss-Lobatto-Runge-Kutta-Methods-1",
    "page": "Runge-Kutta",
    "title": "Gauss-Lobatto Runge-Kutta Methods",
    "category": "section",
    "text": "Function Order Method\ngetTableauLobIIIA2() 2 Gauss-Lobatto IIIA s=2\ngetTableauLobIIIA3() 4 Gauss-Lobatto IIIA s=3\ngetTableauLobIIIA4() 6 Gauss-Lobatto IIIA s=4\ngetTableauLobIIIB2() 2 Gauss-Lobatto IIIB s=2\ngetTableauLobIIIB3() 4 Gauss-Lobatto IIIB s=3\ngetTableauLobIIIB4() 6 Gauss-Lobatto IIIB s=4\ngetTableauLobIIIC2() 2 Gauss-Lobatto IIIC s=2\ngetTableauLobIIIC3() 4 Gauss-Lobatto IIIC s=3\ngetTableauLobIIIC4() 6 Gauss-Lobatto IIIC s=4\ngetTableauLobIIID2() 2 Gauss-Lobatto IIID s=2\ngetTableauLobIIID3() 4 Gauss-Lobatto IIID s=3\ngetTableauLobIIID4() 6 Gauss-Lobatto IIID s=4\ngetTableauLobIIIE2() 2 Gauss-Lobatto IIIE s=2\ngetTableauLobIIIE3() 4 Gauss-Lobatto IIIE s=3\ngetTableauLobIIIE4() 6 Gauss-Lobatto IIIE s=4\ngetTableauLobIIIF2() 4 Gauss-Lobatto IIIF s=2\ngetTableauLobIIIF3() 6 Gauss-Lobatto IIIF s=3\ngetTableauLobIIIF4() 8 Gauss-Lobatto IIIF s=4"
},

{
    "location": "integrators/vprk/#",
    "page": "VPRK",
    "title": "VPRK",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/vprk/#Variational-Partitioned-Runge-Kutta-Integrators-1",
    "page": "VPRK",
    "title": "Variational Partitioned Runge-Kutta Integrators",
    "category": "section",
    "text": ""
},

{
    "location": "integrators/spark/#",
    "page": "SPARK",
    "title": "SPARK",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/spark/#Special-Partitioned-Additive-Runge-Kutta-Integrators-1",
    "page": "SPARK",
    "title": "Special Partitioned Additive Runge-Kutta Integrators",
    "category": "section",
    "text": ""
},

{
    "location": "integrators/cgvi/#",
    "page": "CGVI",
    "title": "CGVI",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/cgvi/#Continuous-Galerkin-Variational-Integrators-1",
    "page": "CGVI",
    "title": "Continuous Galerkin Variational Integrators",
    "category": "section",
    "text": ""
},

{
    "location": "integrators/dgvi/#",
    "page": "DGVI",
    "title": "DGVI",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/dgvi/#Discontinuous-Galerkin-Variational-Integrators-1",
    "page": "DGVI",
    "title": "Discontinuous Galerkin Variational Integrators",
    "category": "section",
    "text": "Discontinuous Galerkin Variational Integrators (DGVIs) are a family of integrators for degenerate Lagrangian systems and for Hamiltonian systems subject to Dirac constraints. For integrators for non-degenerate (regular) Lagrangian and unconstrained Hamiltonian systems see Hamilton-Pontryagin-Galerkin (HPG) Integrators."
},

{
    "location": "integrators/dgvi/#Degenerate-Lagrangian-Systems-1",
    "page": "DGVI",
    "title": "Degenerate Lagrangian Systems",
    "category": "section",
    "text": "Consider a fully degenerate Lagrangian system of the formL(q dotq) = vartheta (q) cdot dotq - H(q) where vartheta (q) denotes the Cartan one-form and H(q) the Hamiltonian, that is usually given by the total energy of the system."
},

{
    "location": "integrators/dgvi/#Discrete-Trajectories-and-Numerical-Quadrature-1",
    "page": "DGVI",
    "title": "Discrete Trajectories and Numerical Quadrature",
    "category": "section",
    "text": "The first step in the derivation of variational integrators is the discretization of the action integralmathcalA q = int limits_0^T L(q(t) dotq(t))  dt To this end, the interval 0T is split into N sub-intervals t_n t_n+1 with t_n = nh and h the time step size, so that t_N = T and the action can be written asmathcalA q = sum limits_n=0^N-1 int limits_t_n t_n+1 L(q(t) dotq(t))  dt Within each interval (t_n t_n+1) a piecewise-polynomial approximation q_h of the trajectory q is constructed using S basis functions varphi_i,q_h(t) vert_(t_n t_n+1) = sum limits_i=1^S x_ni  barvarphi_ni (t) where barvarphi_ni (t) is a rescaled basis function, defined bybarvarphi_ni (t) = varphi_i bigg( fract - t_nt_n+1 - t_n bigg) and it is assumed that varphi_i is compactly supported on 01. These approximations q_h(t) do not need to be continuous across interval boundaries but are indeed allowed to have jumps. Replacing the continuous trajectory q in the action with q_h, we obtainmathcalA q_h = sum limits_n=0^N-1 int limits_(t_n t_n+1) big vartheta (q_h (t)) cdot dotq_h (t) - H(q_h (t)) big  dt\n+ sum limits_n=1^N-1  vartheta (q_h (t)) cdot dotq_h (t) _t=t_n  The integral of the Hamiltonian H(q_h) over the interval boundaries does not contribute to the integral, differently from the term vartheta (q_h) cdot dotq_h, which will determine the numerical flux  cdot _n at t_n of the Discontinuous Galerkin method. The approximation of this term will be discussed below. In order to obtain a fully discrete action, a numerical quadrature rule with R nodes c_i and weights b_i is introduced for the approximation of the integral,mathcalA_d x_d = h sum limits_n=0^N-1 sum limits_i=1^R b_i big vartheta (q_h(t_n + c_i h)) cdot dotq_h (t_n + c_i h) - H(q_h(t_n + c_i h)) big\n+ sum limits_n=1^N-1  vartheta (q_h (t)) cdot dotq_h (t) _t=t_n Here, x_d denotes the vector of all the degrees of freedom, i.e.,x_d = ( x_01  x_0S  x_11  x_N-2S  x_N-11  x_N-1S )^T In order to write the discrete action in a more explicit form, mass and derivative matrices m and a are introduced, whose elements are given bym_ij = varphi_j (c_i) \nqquad\na_ij = varphi_j (c_i) \nqquad\ni = 1  R \n\nj = 1  S so that the solution and its time derivative at the quadrature points can be written asQ_ni equiv q_h(t_n + c_i h) = m_ij x_nj \nqquad\nV_ni equiv dotq_h (t_n + c_i h) = a_ij x_nj wherex_n = ( x_n1  x_nS )^Tis the vector containing the degrees of freedom of q_h vert_t_n t_n+1. Using these definitions, the discrete action can be written asmathcalA_d x_d = h sum limits_n=0^N-1 sum limits_i=1^R b_i big vartheta (Q_ni) cdot V_ni - H(Q_ni) big\n+ sum limits_n=1^N-1  vartheta (q_h (t)) cdot dotq_h (t) _t=t_n "
},

{
    "location": "integrators/dgvi/#Numerical-Fluxes-1",
    "page": "DGVI",
    "title": "Numerical Fluxes",
    "category": "section",
    "text": "In the following, the solution values \"left\" and \"right\" of the jump will be needed. This will be denoted by q_n^- and q_n^+, respectively. Usually, these just correspond to the polynomials on the left and right, evaluated at t_n, i.e.,q_n^- = lim_t uparrow t_n q_h (t) = q_h vert_t_n-1 t_n (t_n) \nqquad\nq_n^+ = lim_t downarrow t_n q_h (t) = q_h vert_t_n t_n+1 (t_n) In principle, however, more general reconstructions of the solution could be used. In the following, it will be assumed, that q_n^- is given by some linear combinations of the degrees of freedom of the polynomial on the left interval and correspondingly that q_n^+ is given by some linear combinations of the degrees of freedom of the polynomial on the right interval, specificallyq_n^- = r^- cdot x_n \nqquad\nq_n^+ = r^+ cdot x_n+1 where r^pm are appropriate coefficient vectors."
},

{
    "location": "integrators/dgvi/#Gauge-Terms-1",
    "page": "DGVI",
    "title": "Gauge Terms",
    "category": "section",
    "text": "The Lagrangian L can be augmented with any total time derivative without changing the (continuous) Euler-Lagrange equations. In particular, one can consider the modified LagrangianL(q dotq) = vartheta (q) cdot dotq - H(q) - nu dfracddt bigg( vartheta (q) cdot q bigg) While this gauge term vanishes in the continuous case, it takes a finite value across jumps of the discontinuous discrete solution, so that the modified discrete action readsmathcalA_d x_d = h sum limits_n=0^N-1 sum limits_i=1^R b_i big vartheta (Q_ni) cdot V_ni - H(Q_ni) big\n+ sum limits_n=1^N-1 biggbigg vartheta (q_h (t)) cdot dotq_h (t) - nu dfracddt bigg( vartheta (q) cdot q bigg) biggbigg_t=t_n In the following, only the modified Lagrangian and action will be considered, in order to obtain a sufficiently general framework for constructing numerical fluxes. For brevity of notation, the prime will be dropped."
},

{
    "location": "integrators/dgvi/#Total-Time-Derivatives-Across-Jumps-1",
    "page": "DGVI",
    "title": "Total Time Derivatives Across Jumps",
    "category": "section",
    "text": "The computation of the total time derivative in the gauge term is simple, at least in the distributional sense. Even though both, vartheta (q_h) and q_h have a jump, the jump occurs at the same position in time, so that the derivative can be computed asdfracddt bigg( vartheta (q_h) cdot q_h bigg) biggvert_t=t_n\n= dfracddt bigg( vartheta (q_n^-) cdot q_n^-  Theta (t_n - t) + vartheta (q_n^+) cdot q_n^+  Theta (t - t_n) bigg) where Theta denotes the Heaviside function. This can be explicitly computed asdfracddt bigg( vartheta (q_h) cdot q_h bigg) biggvert_t=t_n\n= - vartheta (q_n^-) cdot q_n^-  delta (t_n) + vartheta (q_n^+) cdot q_n^+  delta (t_n)\n=  vartheta (q_h) cdot q_h _t=t_n  delta (t_n) with delta (t_n) the Dirac delta-function at t_n."
},

{
    "location": "integrators/dgvi/#Non-conservative-Products-1",
    "page": "DGVI",
    "title": "Non-conservative Products",
    "category": "section",
    "text": "Simple means for integrating the Lagrangian across jumps are provided by discretisations of the integralint limits_0^1 vartheta (Phi(tau q^- q^+)) cdot dfracd Phi(tau q^- q^+)dtau  dtau where Phi is a path connecting the solution values q^- and q^+ on the left and the right of the jump. Upon picking a quadrature rule with sigma nodes gamma_i and corresponding weights beta_i, the discrete product takes the formsum limits_i=1^sigma beta_i  vartheta big( Phi (gamma_i q_n^-  q_n^+) big) cdot dfracdPhidtau (gamma_i q_n^-  q_n^+) For a compact notation, \"mass\" and \"derivative\" vectors mu^pm and alpha^pm are introduced, so thatPhi (gamma_i q_n^-  q_n^+) = mu^-_i q_n^- + mu^+_i q_n^+\nqquad\nPhi (gamma_i q_n^-  q_n^+) = alpha^-_i q_n^- + alpha^+_i q_n^+ and the discrete product can be written assum limits_i=1^sigma beta_i  vartheta ( mu^-_i q_n^- + mu^+_i q_n^+ ) cdot ( alpha^-_i q_n^- + alpha^+_i q_n^+ ) Providing the path Phi by two functions phi^pm(tau), so thatphi(tau q^- q^+) = q^- phi^-(tau) + q^+ phi^+(tau) the components of the \"mass\" and \"derivative\" vectors are given bymu^-_i = phi^- (gamma_i) \nqquad\nmu^+_i = phi^+ (gamma_i) andalpha^-_i = fracdphi^-dtau (gamma_i) \nqquad\nalpha^+_i = fracdphi^+dtau (gamma_i) respectively."
},

{
    "location": "integrators/dgvi/#Discrete-Variational-Principle-1",
    "page": "DGVI",
    "title": "Discrete Variational Principle",
    "category": "section",
    "text": "Using the construction of the previous sections, the discrete action readsmathcalA_d x_d = h sum limits_n=0^N-1 sum limits_i=1^R b_i big vartheta (Q_ni) cdot V_ni - H(Q_ni) big \n+ sum limits_n=1^N-1 bigg sum limits_i=1^sigma beta_i  vartheta ( mu^-_i q_n^- + mu^+_i q_n^+ ) cdot ( alpha^-_i q_n^- + alpha^+_i q_n^+ )\n- nu big vartheta (q_n^+) cdot q_n^+ - vartheta (q_n^-) cdot q_n^- big bigg The discrete Euler-Lagrange equations are obtained by applying Hamilton\'s principle of stationary action to mathcalA_d x_d, that is requiring that delta mathcalA_d x_d = 0. The variations of the discrete action are computed as follows,delta mathcalA_d x_d\n= h sum limits_n=0^N-1 sum limits_i=1^R b_i big delta Q_ni cdot nabla vartheta (Q_ni) cdot V_ni + vartheta (Q_ni) cdot delta V_ni - delta Q_ni cdot nabla H(Q_ni) big \n+ sum limits_n=1^N-1 bigg sum limits_i=1^sigma beta_i  big ( mu^-_i delta q_n^- + mu^+_i delta q_n^+ ) cdot nabla vartheta ( mu^-_i q_n^- + mu^+_i q_n^+ ) cdot ( alpha^-_i q_n^- + alpha^+_i q_n^+ ) + vartheta ( mu^-_i q_n^- + mu^+_i q_n^+ ) cdot ( alpha^-_i delta q_n^- + alpha^+_i delta q_n^+ ) big \n- nu big delta q_n^+ cdot nabla vartheta (q_n^+) cdot q_n^+  + vartheta (q_n^+) cdot delta q_n^+ - delta q_n^- cdot nabla vartheta (q_n^-) cdot q_n^- - vartheta (q_n^-) cdot delta q_n^- big bigg Using the relationsdelta Q_ni = m_ij delta x_nj \nqquad\ndelta V_ni = fraca_ijh delta x_nj \nqquad\ndelta q_n^- = r^-_j delta x_n-1j \nqquad\ndelta q_n^+ = r^+_j delta x_nj the variations of the discrete action becomedelta mathcalA_d x_d\n= sum limits_n=0^N-1 sum limits_i=1^R sum limits_j=1^S b_i big h m_ij delta x_nj cdot nabla vartheta (Q_ni) cdot V_ni + vartheta (Q_ni) cdot a_ij delta x_nj - h m_ij delta x_nj cdot nabla H(Q_ni) big \n+ sum limits_n=1^N-1 bigg sum limits_i=1^sigma sum limits_j=1^S beta_i  big ( mu^-_i r^-_j delta x_n-1j + mu^+_i r^+_j delta x_nj ) cdot nabla vartheta ( mu^-_i q_n^- + mu^+_i q_n^+ ) cdot ( alpha^-_i q_n^- + alpha^+_i q_n^+ ) + vartheta ( mu^-_i q_n^- + mu^+_i q_n^+ ) cdot ( alpha^-_i r^-_j delta x_n-1j + alpha^+_i r^+_j delta x_nj ) big \n- nu big r^+_j delta x_nj cdot nabla vartheta (q_n^+) cdot q_n^+  + vartheta (q_n^+) cdot r^+_j delta x_nj - r^-_j delta x_n-1j cdot nabla vartheta (q_n^-) cdot q_n^- - vartheta (q_n^-) cdot r^-_j delta x_n-1j big bigg Requiring the variation of the discrete action to vanish yields the discrete equations of motion,0 = sum limits_i=1^R b_i big h m_ij nabla vartheta (Q_ni) cdot V_ni + a_ij vartheta (Q_ni) - h m_ij nabla H(Q_ni) big \n+ bigg sum limits_i=1^sigma beta_i  big\n   mu^-_i r^-_j nabla vartheta ( mu^-_i q_n+1^- + mu^+_i q_n+1^+ ) cdot ( alpha^-_i q_n+1^- + alpha^+_i q_n+1^+ )\n + mu^+_i r^+_j nabla vartheta ( mu^-_i q_n  ^- + mu^+_i q_n  ^+ ) cdot ( alpha^-_i q_n  ^- + alpha^+_i q_n  ^+ ) \n + vartheta ( mu^-_i q_n+1^- + mu^+_i q_n+1^+ ) alpha^-_i r^-_j\n + vartheta ( mu^-_i q_n  ^- + mu^+_i q_n  ^+ ) alpha^+_i r^+_j\n   big \n - nu big\n     r^+_j nabla vartheta (q_n  ^+) cdot q_n  ^+ + r^+_j vartheta (q_n  ^+)\n   - r^-_j nabla vartheta (q_n+1^-) cdot q_n+1^- - r^-_j vartheta (q_n+1^-)\n   big\n bigg for all n and all j."
},

{
    "location": "integrators/hpg/#",
    "page": "HPG",
    "title": "HPG",
    "category": "page",
    "text": ""
},

{
    "location": "integrators/hpg/#Hamilton-Pontryagin-Galerkin-Integrators-1",
    "page": "HPG",
    "title": "Hamilton-Pontryagin-Galerkin Integrators",
    "category": "section",
    "text": ""
},

{
    "location": "modules/basis_functions/#",
    "page": "Basis Functions",
    "title": "Basis Functions",
    "category": "page",
    "text": ""
},

{
    "location": "modules/basis_functions/#GeometricIntegrators.BasisFunctions.Basis",
    "page": "Basis Functions",
    "title": "GeometricIntegrators.BasisFunctions.Basis",
    "category": "type",
    "text": "Abstract basis\n\nT: data type   N: number of nodes\n\n\n\n\n\n"
},

{
    "location": "modules/basis_functions/#GeometricIntegrators.BasisFunctions.PolynomialBasis",
    "page": "Basis Functions",
    "title": "GeometricIntegrators.BasisFunctions.PolynomialBasis",
    "category": "type",
    "text": "Abstract polynomial basis.\n\n\n\n\n\n"
},

{
    "location": "modules/basis_functions/#Basis-Functions-1",
    "page": "Basis Functions",
    "title": "Basis Functions",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.BasisFunctions]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/discontinuities/#",
    "page": "Discontinuities",
    "title": "Discontinuities",
    "category": "page",
    "text": ""
},

{
    "location": "modules/discontinuities/#GeometricIntegrators.Discontinuities.PathIntegralLinear",
    "page": "Discontinuities",
    "title": "GeometricIntegrators.Discontinuities.PathIntegralLinear",
    "category": "type",
    "text": "PathIntegralLinear is a path integral along a linear path\n\nphi (tau q^- q^+) = (1-tau) q^- + tau q^+ \n\n\n\n\n\n"
},

{
    "location": "modules/discontinuities/#GeometricIntegrators.Discontinuities.PathIntegralTrigonometric",
    "page": "Discontinuities",
    "title": "GeometricIntegrators.Discontinuities.PathIntegralTrigonometric",
    "category": "type",
    "text": "PathIntegralTrigonometric is a path integral along a cos^2/sin^2 path\n\nphi (tau q^- q^+) = cos^2 (pi tau  2) q^- + sin^2 (pi tau  2) q^+ \n\n\n\n\n\n"
},

{
    "location": "modules/discontinuities/#Discontinuities-1",
    "page": "Discontinuities",
    "title": "Discontinuities",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Discontinuities]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/equations/#",
    "page": "Equations",
    "title": "Equations",
    "category": "page",
    "text": ""
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.DAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.DAE",
    "category": "type",
    "text": "DAE: Differential Algebraic Equation\n\nDefines a differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t)) + u(t q(t) lambda(t))   q(t_0) = q_0  \n0 = phi (t q(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector field v, projection u, algebraic constraint phi=0, initial conditions q_0 and lambda_0, the dynamical variable q taking values in mathbbR^m and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nd: dimension of dynamical variable q and the vector field v\nm: dimension of algebraic variable lambda and the constraint phi\nn: number of initial conditions\nv: function computing the vector field\nu: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\nλ₀: initial condition for algebraic variable lambda\n\nThe function v, providing the vector field, takes three arguments, v(t, q, v), the functions u and ϕ, providing the projection and the algebraic constraint take four arguments, u(t, q, λ, u) and ϕ(t, q, λ, ϕ), where t is the current time, q and λ are the current solution vectors, and v, u and ϕ are the vectors which hold the result of evaluating the vector field v, the projection u and the algebraic constraint phi on t, q and λ.\n\nExample\n\n    function v(t, q, v)\n        v[1] = q[1]\n        v[2] = q[2]\n    end\n\n    function u(t, q, λ, u)\n        u[1] = +λ[1]\n        u[2] = -λ[1]\n    end\n\n    function ϕ(t, q, λ, ϕ)\n        ϕ[1] = q[2] - q[1]\n    end\n\n    t₀ = 0.\n    q₀ = [1., 1.]\n    λ₀ = [0.]\n\n    dae = DAE(v, u, ϕ, t₀, q₀, λ₀)\n\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.HDAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.HDAE",
    "category": "type",
    "text": "HDAE: Hamiltonian Differential Algebraic Equation\n\nDefines a Hamiltonian differential algebraic initial value problem, that is a canonical Hamiltonian system of equations subject to Dirac constraints,\n\nbeginalign*\ndotq (t) = v_1(t q(t) p(t)) + v_2(t q(t) p(t) lambda(t)) + v_3(t q(t) p(t) lambda(t) gamma(t))   q(t_0) = q_0  \ndotp (t) = f_1(t q(t) p(t)) + f_2(t q(t) p(t) lambda(t)) + f_3(t q(t) p(t) lambda(t) gamma(t))   p(t_0) = p_0  \n0 = phi (t q(t) p(t))  \n0 = psi (t q(t) p(t) lambda(t)) \nendalign*\n\nwith vector fields v_i and f_i for i = 1  3, primary constraint phi(qp)=0 and secondary constraint psi(qplambda)=0, initial conditions (q_0 p_0), the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d and the algebraic variables (lambda gamma) taking values in mathbbR^n times mathbbR^d.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields v and f\nm: dimension of algebraic variables lambda and gamma and the constraints phi and psi\nn: number of initial conditions\nv: tuple of functions computing the vector fields v_i, i = 1  3\nf: tuple of functions computing the vector fields f_i, i = 1  3\nϕ: primary constraints\nψ: secondary constraints\nt₀: initial time\nq₀: initial condition for dynamical variable q\np₀: initial condition for dynamical variable p\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.IDAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.IDAE",
    "category": "type",
    "text": "IDAE: Implicit Differential Algebraic Equation\n\nDefines a partitioned differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = v(t) + u(t q(t) p(t) lambda(t))   q(t_0) = q_0  \ndotp (t) = f(t q(t) v(t)) + r(t q(t) p(t) lambda(t))   p(t_0) = p_0  \np(t) = p(t q(t) v(t))   \n0 = phi (t q(t) p(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector field f, the momentum defined by p, projection u and r, algebraic constraint phi=0, conditions (q_0 p_0) and lambda_0, the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and p\nm: dimension of algebraic variable lambda and the constraint phi\nn: number of initial conditions\nf: function computing the vector field f\np: function computing p\nu: function computing the projection\ng: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\np₀: initial condition for dynamical variable p\nλ₀: initial condition for algebraic variable lambda\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.IODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.IODE",
    "category": "type",
    "text": "IODE: Implicit Ordinary Differential Equation\n\nDefines an implicit initial value problem\n\nbeginalign*\ndotq (t) = v(t)  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) v(t))  \np(t_0) = p_0  \np(t) = α(t q(t) v(t))\nendalign*\n\nwith vector field f, the momentum defined by p, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d. This is a special case of a differential algebraic equation with dynamical variables (qp) and algebraic variable v.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and p\nα: function determining the momentum\nf: function computing the vector field\ng: function determining the projection, given by ∇α(q)λ\nv: function computing an initial guess for the velocity field (optional)\nt₀: initial time (optional)\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions α and f must have the interface\n\n    function α(t, q, v, p)\n        p[1] = ...\n        p[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, v, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, v is the current velocity and f and p are the vectors which hold the result of evaluating the functions f and α on t, q and v. The funtions g and v are specified by\n\n    function g(t, q, λ, g)\n        g[1] = ...\n        g[2] = ...\n        ...\n    end\n\nand\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.ODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.ODE",
    "category": "type",
    "text": "ODE: Ordinary Differential Equation\n\nDefines an initial value problem\n\ndotq (t) = v(t q(t))  qquad q(t_0) = q_0 \n\nwith vector field v, initial condition q_0 and the solution q taking values in mathbbR^d.\n\nFields\n\nd: dimension of dynamical variable q and the vector field v\nn: number of initial conditions\nv: function computing the vector field\nt₀: initial time\nq₀: initial condition\n\nThe function v providing the vector field must have the interface\n\n    function v(t, q, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, and v is the vector which holds the result of evaluating the vector field v on t and q.\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.PDAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PDAE",
    "category": "type",
    "text": "PDAE: Partitioned Differential Algebraic Equation\n\nDefines a partitioned differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t) p(t)) + u(t q(t) p(t) lambda(t))   q(t_0) = q_0  \ndotp (t) = f(t q(t) p(t)) + r(t q(t) p(t) lambda(t))   p(t_0) = p_0  \n0 = phi (t q(t) p(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector fields v and f, projection u and r, algebraic constraint phi=0, conditions (q_0 p_0) and lambda_0, the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and p\nm: dimension of algebraic variable lambda and the constraint phi\nn: number of initial conditions\nv: function computing the vector field v\nf: function computing the vector field f\nu: function computing the projection\ng: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\np₀: initial condition for dynamical variable p\nλ₀: initial condition for algebraic variable lambda\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.PODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PODE",
    "category": "type",
    "text": "IODE: Partitioned Ordinary Differential Equation\n\nDefines a partitioned initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t) p(t))  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) p(t))  \np(t_0) = p_0 \nendalign*\n\nwith vector fields v and f, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields v and f\nv: function computing the vector field v\nf: function computing the vector field f\nt₀: initial time\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions v and f must have the interface\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, p, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q and p are the current solution vectors and v and f are the vectors which hold the result of evaluating the vector fields v and f on t, q and p.\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.PSDE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PSDE",
    "category": "type",
    "text": "PSDE: Stratonovich Partitioned Stochastic Differential Equation\n\nDefines a partitioned stochastic differential initial value problem\n\nbeginalign*\ndq (t) = v(t q(t))  dt + B(t q(t)) circ dW   q(t_0) = q_0 \ndp (t) = f(t q(t))  dt + G(t q(t)) circ dW   p(t_0) = p_0\nendalign*\n\nwith the drift vector fields v and f, diffusion matrices B and G, initial conditions q_0 and p_0, the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d, and the m-dimensional Wiener process W\n\nFields\n\nd:  dimension of dynamical variable q and the vector field v\nm:  dimension of the Wiener process\nn:  number of initial conditions\nns: number of sample paths\nv:  function computing the drift vector field for the position variable q\nf:  function computing the drift vector field for the momentum variable p\nB:  function computing the d x m diffusion matrix for the position variable q\nG:  function computing the d x m diffusion matrix for the momentum variable p\nt₀: initial time\nq₀: initial condition for dynamical variable q (may be a random variable itself)\np₀: initial condition for dynamical variable p (may be a random variable itself)\n\nThe functions v, f, \'B\' and G, providing the drift vector fields and diffusion matrices, take four arguments, v(t, q, p, v), f(t, q, p, f), B(t, q, p,  B) and G(t, q, p, G), where t is the current time, (q, p) is the current solution vector, and v, f, \'B\' and G are the variables which hold the result of evaluating the vector fields v, f and the matrices B, G on t and (q,p).\n\nExample\n\n    function v(λ, t, q, v)\n        v[1] = λ*q[1]\n        v[2] = λ*q[2]\n    end\n\n    function B(μ, t, q, B)\n        B[1] = μ*q[1]\n        B[2] = μ*q[2]\n    end\n\n    t₀ = 0.\n    q₀ = [1., 1.]\n    λ  = 2.\n    μ  = 1.\n\n    v_sde = (t, q, v) -> v(λ, t, q, v)\n    B_sde = (t, q, B) -> B(μ, t, q, B)\n\n    sde = SDE(v_sde, B_sde, t₀, q₀)\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.SDE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.SDE",
    "category": "type",
    "text": "SDE: Stratonovich Stochastic Differential Equation\n\nDefines a stochastic differential initial value problem\n\nbeginalign*\ndq (t) = v(t q(t))  dt + B(t q(t)) circ dW   q(t_0) = q_0 \nendalign*\n\nwith drift vector field v, diffusion matrix B, initial conditions q_0, the dynamical variable q taking values in mathbbR^d, and the m-dimensional Wiener process W\n\nFields\n\nd:  dimension of dynamical variable q and the vector field v\nm:  dimension of the Wiener process\nn:  number of initial conditions\nns: number of sample paths\nv:  function computing the deterministic vector field\nB:  function computing the d x m diffusion matrix\nt₀: initial time\nq₀: initial condition for dynamical variable q (may be a random variable itself)\n\nThe functions v and B, providing the drift vector field and diffusion matrix, v(t, q, v) and B(t, q, B; col=0), where t is the current time, q is the current solution vector, and v and B are the variables which hold the result of evaluating the vector field v and the matrix B on t and q (if col==0), or the column col of the matrix B (if col>0).\n\nExample\n\n    function v(λ, t, q, v)\n        v[1] = λ*q[1]\n        v[2] = λ*q[2]\n    end\n\n    function B(μ, t, q, B; col=0)\n        if col==0 #whole matrix\n            B[1,1] = μ*q[1]\n            B[2,1] = μ*q[2]\n        elseif col==1\n            #just first column\n        end\n    end\n\n    t₀ = 0.\n    q₀ = [1., 1.]\n    λ  = 2.\n    μ  = 1.\n\n    v_sde = (t, q, v) -> v(λ, t, q, v)\n    B_sde = (t, q, B) -> B(μ, t, q, B)\n\n    sde = SDE(v_sde, B_sde, t₀, q₀)\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.SODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.SODE",
    "category": "type",
    "text": "SODE: Split Ordinary Differential Equation\n\nDefines an initial value problem\n\ndotq (t) = v(t q(t))  qquad q(t_0) = q_0 \n\nwith vector field v, initial condition q_0 and the solution q taking values in mathbbR^d. Here, the vector field v is given as a sum of vector fields\n\nv (t) = v_1 (t) +  + v_r (t) \n\nFields\n\nd: dimension of dynamical variable q and the vector field v\nv: tuple of functions computing the vector field\nt₀: initial time\nq₀: initial condition\n\nThe functions v_i providing the vector field must have the interface\n\n    function v_i(t, q₀, q₁, h)\n        q₁[1] = q₀[1] + ...\n        q₁[2] = q₀[2] + ...\n        ...\n    end\n\nwhere t is the current time, q₀ is the current solution vector, q₁ is the new solution vector which holds the result of computing one substep with the vector field v_i on t and q₀, and h is the (sub-)timestep to compute the update for.\n\nThe fact that the function v returns the solution and not just the vector field for each substep increases the flexibility for the use of splitting methods, e.g., it allows to use another integrator for solving substeps.\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.SPSDE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.SPSDE",
    "category": "type",
    "text": "SPSDE: Stratonovich Split Partitioned Stochastic Differential Equation\n\nDefines a partitioned stochastic differential initial value problem\n\nbeginalign*\ndq (t) =   v(t q(t))  dt + B(t q(t)) circ dW   q(t_0) = q_0 \ndp (t) =  f1(t q(t)) + f2(t q(t))   dt +  G1(t q(t)) + G2(t q(t))  circ dW   p(t_0) = p_0\nendalign*\n\nwith the drift vector fields v and fi, diffusion matrices B and Gi, initial conditions q_0 and p_0, the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d, and the m-dimensional Wiener process W\n\nFields\n\nd:  dimension of dynamical variable q and the vector fields vi\nm:  dimension of the Wiener process\nn:  number of initial conditions\nns: number of sample paths\nv :  function computing the drift vector field for the position variable q\nf1:  function computing the drift vector field for the momentum variable p\nf2:  function computing the drift vector field for the momentum variable p\nB :  function computing the d x m diffusion matrix for the position variable q\nG1:  function computing the d x m diffusion matrix for the momentum variable p\nG2:  function computing the d x m diffusion matrix for the momentum variable p\nt₀: initial time\nq₀: initial condition for dynamical variable q (may be a random variable itself)\np₀: initial condition for dynamical variable p (may be a random variable itself)\n\nThe functions v, f, \'B\' and G, providing the drift vector fields and diffusion matrices, take four arguments, v(t, q, p, v), f(t, q, p, f), B(t, q, p,  B) and G(t, q, p, G), where t is the current time, (q, p) is the current solution vector, and v, f, \'B\' and G are the variables which hold the result of evaluating the vector fields v, f and the matrices B, G on t and (q,p).\n\nExample\n\n    function v(λ, t, q, v)\n        v[1] = λ*q[1]\n        v[2] = λ*q[2]\n    end\n\n    function B(μ, t, q, B)\n        B[1] = μ*q[1]\n        B[2] = μ*q[2]\n    end\n\n    t₀ = 0.\n    q₀ = [1., 1.]\n    λ  = 2.\n    μ  = 1.\n\n    v_sde = (t, q, v) -> v(λ, t, q, v)\n    B_sde = (t, q, B) -> B(μ, t, q, B)\n\n    sde = SDE(v_sde, B_sde, t₀, q₀)\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#GeometricIntegrators.Equations.VODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.VODE",
    "category": "type",
    "text": "VODE: Variational Ordinary Differential Equation\n\nDefines an implicit initial value problem\n\nbeginalign*\ndotq (t) = v(t)  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) v(t))  \np(t_0) = p_0  \np(t) = α(t q(t) v(t))\nendalign*\n\nwith vector field f, the momentum defined by p, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d. This is a special case of a differential algebraic equation with dynamical variables (qp) and algebraic variable v.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and p\nα: function determining the momentum\nf: function computing the vector field\ng: function determining the projection, given by ∇α(q)λ\nv: function computing an initial guess for the velocity field (optional)\nt₀: initial time (optional)\nq₀: initial condition for q\np₀: initial condition for p\nλ₀: initial condition for λ\n\nThe functions α and f must have the interface\n\n    function α(t, q, v, p)\n        p[1] = ...\n        p[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, v, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, v is the current velocity and f and p are the vectors which hold the result of evaluating the functions f and α on t, q and v. The funtions g and v are specified by\n\n    function g(t, q, λ, g)\n        g[1] = ...\n        g[2] = ...\n        ...\n    end\n\nand\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\n\n\n\n\n"
},

{
    "location": "modules/equations/#Equations-1",
    "page": "Equations",
    "title": "Equations",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Equations]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/integrators/#",
    "page": "Integrators",
    "title": "Integrators",
    "category": "page",
    "text": ""
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.AbstractTableau",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.AbstractTableau",
    "category": "type",
    "text": "Holds the information for the various methods\' tableaus.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.AbstractTableauIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.AbstractTableauIRK",
    "category": "type",
    "text": "Holds the tableau of an implicit Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.AbstractTableauPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.AbstractTableauPRK",
    "category": "type",
    "text": "Holds the tableau of a partitioned Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.AbstractTableauRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.AbstractTableauRK",
    "category": "type",
    "text": "Holds the tableau of a Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.CoefficientsARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.CoefficientsARK",
    "category": "type",
    "text": "Holds the coefficients of an additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.CoefficientsMRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.CoefficientsMRK",
    "category": "type",
    "text": "Holds the multiplier Runge-Kutta coefficients.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.CoefficientsPGLRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.CoefficientsPGLRK",
    "category": "type",
    "text": "Holds the coefficients of a projected Gauss-Legendre Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.CoefficientsPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.CoefficientsPRK",
    "category": "type",
    "text": "Holds the coefficients of a projective Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.CoefficientsRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.CoefficientsRK",
    "category": "type",
    "text": "Holds the coefficients of a Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{DAE,TableauARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for additive Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{DAE,TableauSARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for special additive Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{Equation,AbstractTableau,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Print error for integrators not implemented, yet.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{IDAE,TableauVPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for variational partitioned additive Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{IDAE,TableauVSPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for variational special partitioned additive Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{IODE,CoefficientsPGLRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for Projected Gauss-Legendre Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{IODE,TableauVPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for variational partitioned Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{ODE,TableauDIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for diagonally implicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{ODE,TableauERK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for explicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{ODE,TableauFIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for fully implicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{PDAE,TableauPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for partitioned additive Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{PDAE,TableauSPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for special partitioned additive Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{PODE,TableauEPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for explicit partitioned Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{PODE,TableauIPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for implicit partitioned Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{PSDE,TableauSIPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for stochastic fully implicit partitioned Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{SDE,TableauSERK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for stochastic explicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{SDE,TableauSIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for stochastic fully implicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{SDE,TableauWERK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for weak explicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{SDE,TableauWIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for weak fully implicit Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{SODE,AbstractTableauSplitting,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for splitting tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{SPSDE,TableauSISPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for stochastic fully implicit split partitioned Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.Integrator-Tuple{VODE,TableauFIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "method",
    "text": "Create integrator for formal Lagrangian Runge-Kutta tableau.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorCGVI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorCGVI",
    "category": "type",
    "text": "Continuous Galerkin Variational Integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorDGVI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDGVI",
    "category": "type",
    "text": "IntegratorDGVI: Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nFields\n\nequation: Implicit Ordinary Differential Equation\nbasis: piecewise polynomial basis\nquadrature: numerical quadrature rule\nΔt: time step\nparams: ParametersDGVI\nsolver: nonlinear solver\niguess: initial guess\nq: current solution vector for trajectory\np: current solution vector for one-form\ncache: temporary variables for nonlinear solver\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorDGVIEXP",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDGVIEXP",
    "category": "type",
    "text": "IntegratorDGVIEXP: Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nFields\n\nequation: Implicit Ordinary Differential Equation\nbasis: piecewise polynomial basis\nquadrature: numerical quadrature rule\nΔt: time step\nparams: ParametersDGVIEXP\nsolver: nonlinear solver\niguess: initial guess\nq: current solution vector for trajectory\np: current solution vector for one-form\ncache: temporary variables for nonlinear solver\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorDGVIP0",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDGVIP0",
    "category": "type",
    "text": "IntegratorDGVIP0: Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nFields\n\nequation: Implicit Ordinary Differential Equation\nbasis: piecewise polynomial basis\nquadrature: numerical quadrature rule\nΔt: time step\nparams: ParametersDGVIP0\nsolver: nonlinear solver\niguess: initial guess\nq: current solution vector for trajectory\np: current solution vector for one-form\ncache: temporary variables for nonlinear solver\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorDGVIP1",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDGVIP1",
    "category": "type",
    "text": "IntegratorDGVIP1: Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nFields\n\nequation: Implicit Ordinary Differential Equation\nbasis: piecewise polynomial basis\nquadrature: numerical quadrature rule\nΔt: time step\nparams: ParametersDGVIP1\nsolver: nonlinear solver\niguess: initial guess\nq: current solution vector for trajectory\np: current solution vector for one-form\ncache: temporary variables for nonlinear solver\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorDGVIPI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDGVIPI",
    "category": "type",
    "text": "IntegratorDGVIPI: Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nFields\n\nequation: Implicit Ordinary Differential Equation\nbasis: piecewise polynomial basis\nquadrature: numerical quadrature rule\njump: jump across discontinuity\nΔt: time step\nparams: ParametersDGVIPI\nsolver: nonlinear solver\niguess: initial guess\nq: current solution vector for trajectory\nq⁻: current solution vector for trajectory, lhs of jump\nq⁺: current solution vector for trajectory, rhs of jump\ncache: temporary variables for nonlinear solver\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorDIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDIRK",
    "category": "type",
    "text": "Diagonally implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorEPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorEPRK",
    "category": "type",
    "text": "Explicit partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorERK",
    "category": "type",
    "text": "Explicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorFIRK",
    "category": "type",
    "text": "Fully implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorFLRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorFLRK",
    "category": "type",
    "text": "Formal Lagrangian Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorGPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorGPARK",
    "category": "type",
    "text": "Special Partitioned Additive Runge Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorHSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorHSPARK",
    "category": "type",
    "text": "Hamiltonian Specialised Partitioned Additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorIPRK",
    "category": "type",
    "text": "Implicit partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorPARK",
    "category": "type",
    "text": "Implicit partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorPGLRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorPGLRK",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSERK",
    "category": "type",
    "text": "Stochastic Explicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSIPRK",
    "category": "type",
    "text": "Stochastic implicit partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSIRK",
    "category": "type",
    "text": "Stochastic implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSISPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSISPRK",
    "category": "type",
    "text": "Stochastic implicit partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSPARK",
    "category": "type",
    "text": "Variational special partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSplitting",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSplitting",
    "category": "type",
    "text": "Splitting integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSplitting-Union{Tuple{ST}, Tuple{VT}, Tuple{TT}, Tuple{DT}, Tuple{SODE{DT,TT,VT,N} where N,ST,TT}} where ST<:TableauSplittingGS{TT} where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSplitting",
    "category": "method",
    "text": "Construct splitting integrator for symmetric splitting tableau with general stages.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSplitting-Union{Tuple{ST}, Tuple{VT}, Tuple{TT}, Tuple{DT}, Tuple{SODE{DT,TT,VT,N} where N,ST,TT}} where ST<:TableauSplittingNS{TT} where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSplitting",
    "category": "method",
    "text": "Construct splitting integrator for non-symmetric splitting tableau with general stages.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorSplitting-Union{Tuple{ST}, Tuple{VT}, Tuple{TT}, Tuple{DT}, Tuple{SODE{DT,TT,VT,N} where N,ST,TT}} where ST<:TableauSplittingSS{TT} where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSplitting",
    "category": "method",
    "text": "Construct splitting integrator for symmetric splitting tableau with symmetric stages.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPARK",
    "category": "type",
    "text": "Variational partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRK",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKdegenerate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKdegenerate",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKpLegendre",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpLegendre",
    "category": "type",
    "text": "Variational special partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKpMidpoint",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpMidpoint",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKpSecondary",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpSecondary",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator with projection on secondary constraint.\n\nbeginalign*\nP_ni = dfracpartial Lpartial v (Q_ni V_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  big( V_nj + Lambda_nj big)  \nq_n+1 = q_n + h sum limits_i=1^s b_i  big( V_ni + Lambda_ni big)  \nF_ni = dfracpartial Lpartial q (Q_ni V_ni)  \nP_ni = p_n + h sum limits_i=1^s bara_ij  big( F_nj + nabla vartheta (Q_nj) cdot Lambda_nj big) - d_i lambda  \np_n+1 = p_n + h sum limits_i=1^s barb_i  big( F_ni + nabla vartheta (Q_nj) cdot Lambda_nj big)  \n0 = sum limits_i=1^s d_i V_i  \n0 = sum limits_j=1^s omega_ij Psi_nj  \n0 = phi (q_n+1 p_n+1) \nendalign*\n\nsatisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + barb_j a_ji = b_i barb_j  \nbarb_i = b_i \nendalign*\n\nthe primary constraint,\n\nbeginalign*\nphi(qp) = p - vartheta (q) = 0 \nendalign*\n\nat the final solution (q_n+1 p_n+1), and super positions of the secondary constraints,\n\nbeginalign*\npsi(qdotqpdotp)\n= dotp - dotq cdot nabla vartheta (q)\n= big( nabla vartheta (q) - nabla vartheta^T (q) big) cdot dotq - nabla H (q)\n= 0\nendalign*\n\nat the internal stages,\n\nbeginalign*\nPsi_nj = big( nabla vartheta (Q_nj) - nabla vartheta^T (Q_nj) big) cdot V_nj - nabla H (Q_nj) \nendalign*\n\nHere, omega is a (s-1) times s matrix, chosen such that the resulting method has optimal order. The vector d is zero for Gauss-Legendre methods and needs to be chosen appropriately for Gauss-Lobatto methods (for details see documentation of VPRK methods).\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKpStandard",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpStandard",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKpSymmetric",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpSymmetric",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVPRKpVariational",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpVariational",
    "category": "type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorVSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVSPARK",
    "category": "type",
    "text": "Variational Specialised Partitioned Additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorWERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorWERK",
    "category": "type",
    "text": "Stochastic Explicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.IntegratorWIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorWIRK",
    "category": "type",
    "text": "Stochastic implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauARK",
    "category": "type",
    "text": "Holds the tableau of a additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauDIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauDIRK",
    "category": "type",
    "text": "Holds the tableau of a diagonally implicit Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauEPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauEPRK",
    "category": "type",
    "text": "TableauEPRK: Tableau of an Explicit Partitioned Runge-Kutta method\n\nbeginalign*\nV_ni = hphantom- dfracpartial Hpartial p (Q_ni P_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  V_nj  \nq_n+1 = q_n + h sum limits_i=1^s b_i  V_ni  \nF_ni = - dfracpartial Hpartial q (Q_ni P_ni)  \nP_ni = p_n + h  sum limits_i=1^s bara_ij  F_nj  \np_n+1 = p_n + h sum limits_i=1^s barb_i  F_ni \nendalign*\n\nusually satisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + b_j a_ji = b_i b_j  \nbarb_i = b_i \nendalign*\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauERK",
    "category": "type",
    "text": "Holds the tableau of an explicit Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauFIRK",
    "category": "type",
    "text": "Holds the tableau of a fully implicit Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauGLM",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauGLM",
    "category": "type",
    "text": "Holds the tableau of a general linear method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauGPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauGPARK",
    "category": "type",
    "text": "Holds the tableau of a spezialized partitioned additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauHSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauHSPARK",
    "category": "type",
    "text": "Holds the tableau of an Hamiltonian Specialised Partitioned Additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauIPRK",
    "category": "type",
    "text": "TableauIPRK: Tableau of an Implicit Partitioned Runge-Kutta method\n\nbeginalign*\nP_ni = dfracpartial Lpartial v (Q_ni V_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  V_nj  \nq_n+1 = q_n + h sum limits_i=1^s b_i  V_ni  \nF_ni = dfracpartial Lpartial q (Q_ni V_ni)  \nP_ni = p_n + h  sum limits_i=1^s bara_ij  F_nj  \np_n+1 = p_n + h sum limits_i=1^s barb_i  F_ni \nendalign*\n\nusually satisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + b_j a_ji = b_i b_j  \nbarb_i = b_i \nendalign*\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauPARK",
    "category": "type",
    "text": "Holds the tableau of an partitioned additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSARK",
    "category": "type",
    "text": "Holds the tableau of a spezialized additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSERK",
    "category": "type",
    "text": "Holds the tableau of a stochastic explicit Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSIPRK",
    "category": "type",
    "text": "Holds the tableau of a stochastic implicit partitioned Runge-Kutta method.  qdrift, pdrift hold the RK coefficients for the drift part,  and qdiff, pdiff hold the RK coefficients for the diffusion part of the SDE.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSIRK",
    "category": "type",
    "text": "Holds the tableau of a stochastic implicit Runge-Kutta method.  qdrift holds the RK coefficients for the drift part,  and qdiff holds the RK coefficients for the diffusion part of the SDE.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSISPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSISPRK",
    "category": "type",
    "text": "Holds the tableau of a stochastic implicit split partitioned Runge-Kutta method.  qdrift, pdrift1, pdrift2 hold the RK coefficients for the drift parts,  and qdiff, pdiff1, pdiff2 hold the RK coefficients for the diffusion part of the SDE.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSPARK",
    "category": "type",
    "text": "Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSplittingGS",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSplittingGS",
    "category": "type",
    "text": "Tableau for symmetric splitting methods with general stages.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSplittingNS",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSplittingNS",
    "category": "type",
    "text": "Tableau for non-symmetric splitting methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauSplittingSS",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSplittingSS",
    "category": "type",
    "text": "Tableau for symmetric splitting methods with symmetric stages.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauVPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauVPARK",
    "category": "type",
    "text": "Holds the tableau of an variational partitioned additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauVPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauVPRK",
    "category": "type",
    "text": "TableauVPRK: Tableau of a Variational Partitioned Runge-Kutta method\n\nbeginalign*\nP_ni = dfracpartial Lpartial v (Q_ni V_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  V_nj  \nq_n+1 = q_n + h sum limits_i=1^s b_i  V_ni  \nF_ni = dfracpartial Lpartial q (Q_ni V_ni)  \nP_ni = p_n + h sum limits_i=1^s bara_ij  F_nj - d_i lambda  \np_n+1 = p_n + h sum limits_i=1^s barb_i  F_ni  \n\n0 = sum limits_i=1^s d_i V_i  \nendalign*\n\nsatisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + b_j a_ji = b_i b_j  \nbarb_i = b_i \nendalign*\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauVSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauVSPARK",
    "category": "type",
    "text": "Holds the tableau of an Variational Specialised Partitioned Additive Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauWERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauWERK",
    "category": "type",
    "text": "Holds the tableau of a weak explicit Runge-Kutta method.\n\nAccording to Andreas Rossler, \"Second order Runge-Kutta methods for Stratonovich stochastic differential equations\",    BIT Numerical Mathematics (2007) 47, equation (5.1)\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.TableauWIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauWIRK",
    "category": "type",
    "text": "Holds the tableau of a weak implicit Runge-Kutta method.\n\nAccording to Wang, Hong, Xu, \"Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems\",    Commun. Comput. Phys. 21(1), 2017\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{FR}, Tuple{QR}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersDGVIPI{DT,TT,D,S,QR,FR,ΘT,FT,GT} where GT where FT where ΘT}} where FR where QR where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersCGVI{DT,TT,D,S,R,ΘT,FT} where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersDGVIEXP{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersDGVIP0{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersDGVIP1{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersDGVI{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersDIRK{DT,TT,ET,D,S},Int64}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of fully implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersFIRK{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of fully implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersIPRK{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of implicit partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKdegenerate{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute solution of degenerate symplectic partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpInternal{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpMidpoint{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpSecondary{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpStandard{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of projected variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpSymmetric{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpVariational{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of projected variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRK{DT,TT,ET,D,S}}} where S where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{FT}, Tuple{ΘT}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersVPRKpLegendre{DT,TT,ΘT,FT,D,S}}} where S where D where FT where ΘT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational special partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersPGLRK{DT,TT,D,S,ET} where ET}} where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{D}, Tuple{VT}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersFLRK{DT,TT,VT,D,S}}} where S where D where VT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of formal Lagrangian Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersSIPRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersSIRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of stochastic implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersSISPRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of stochastic implicit split partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{ST,1},ParametersWIRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of weak implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{ψT}, Tuple{ϕT}, Tuple{FT}, Tuple{VT}, Tuple{TT}, Tuple{DT}, Tuple{Array{DT,1},Array{DT,1},ParametersHSPARK{DT,TT,VT,FT,ϕT,ψT}}} where ψT where ϕT where FT where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of Hamiltonian Specialised Partitioned Additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{ϕT}, Tuple{GT}, Tuple{UT}, Tuple{PT}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{Array{DT,1},Array{DT,1},ParametersPARK{DT,TT,FT,PT,UT,GT,ϕT}}} where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{ϕT}, Tuple{GT}, Tuple{UT}, Tuple{PT}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{Array{DT,1},Array{DT,1},ParametersSPARK{DT,TT,FT,PT,UT,GT,ϕT}}} where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational special partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.function_stages!-Union{Tuple{ϕT}, Tuple{GT}, Tuple{UT}, Tuple{PT}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{Array{DT,1},Array{DT,1},ParametersVPARK{DT,TT,FT,PT,UT,GT,ϕT}}} where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.function_stages!",
    "category": "method",
    "text": "Compute stages of variational partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initialize!-Tuple{GeometricIntegrators.Integrators.StochasticIntegrator,StochasticSolution,Int64,Int64,Int64,Int64}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initialize!",
    "category": "method",
    "text": "Initialize stochastic integrator for the sample paths k with k₁ ≤ k ≤ k₂, initial conditions m with m₁ ≤ m ≤ m₂ and time step 0.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initialize!-Tuple{Integrator,Solution,Int64,Int64}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initialize!",
    "category": "method",
    "text": "Initialize integrator for initial conditions m with m₁ ≤ m ≤ m₂ and time step 0.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "function",
    "text": "Apply integrator for ntime time steps and return solution.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "function",
    "text": "Integrate given equation with given tableau for ntime time steps and return solution.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "function",
    "text": "Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.StochasticIntegrator,StochasticSolution,Int64,Int64,Int64,Int64}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate SDE for the sample paths k with k₁ ≤ k ≤ k₂ and initial conditions m with m₁ ≤ m ≤ m₂.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.StochasticIntegrator,StochasticSolution}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate SDE for all sample paths and initial conditions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Tuple{Integrator,Solution,Any,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Tuple{Integrator,Solution}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate ODE for all initial conditions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Tuple{IntegratorGPARK,SolutionPDAE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate partitioned DAE with Special Additive Runge Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Union{Tuple{N}, Tuple{TT}, Tuple{DT}, Tuple{Integrator{DT,TT},Solution{DT,TT,N},Int64,Int64,Int64,Int64}} where N where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂ for time steps n with n₁ ≤ n ≤ n₂.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate!-Union{Tuple{N}, Tuple{TT}, Tuple{DT}, Tuple{StochasticIntegrator{DT,TT},StochasticSolution{DT,TT,N,NW} where NW,Int64,Int64,Int64,Int64,Int64,Int64}} where N where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "method",
    "text": "Integrate SDE for the sample paths k with k₁ ≤ k ≤ k₂, the initial conditions m with m₁ ≤ m ≤ m₂, and the time steps n with n₁ ≤ n ≤ n₂.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.readTableauERKFromFile-Tuple{AbstractString,AbstractString}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.readTableauERKFromFile",
    "category": "method",
    "text": "Read explicit Runge-Kutta tableau from file.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.writeTableauToFile-Union{Tuple{T}, Tuple{AbstractString,AbstractTableauRK{T}}} where T",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.writeTableauToFile",
    "category": "method",
    "text": "Write Runge-Kutta tableau to file.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.AbstractTableauERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.AbstractTableauERK",
    "category": "type",
    "text": "Holds the tableau of an explicit Runge-Kutta method.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.NonlinearFunctionCacheDGVI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionCacheDGVI",
    "category": "type",
    "text": "Nonlinear function cache for Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nST: data type\nD: number of dimensions\nS: number of degrees of freedom\nR: number of nodes of quadrature formula\n\nFields\n\nX: degrees of freedom\nQ: solution at quadrature nodes\nV: velocity at quadrature nodes\nP: one-form at quadrature nodes\nF: forces at quadrature nodes\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nq⁺: current solution of qₙ⁺\nq̅:  current solution of qₙ₊₁\nq̅⁻: current solution of qₙ₊₁⁻\nq̅⁺: current solution of qₙ₊₁⁺\nϕ:  average of the solution at tₙ\nϕ̅:  average of the solution at tₙ₊₁\nλ:  jump of the solution at tₙ\nλ̅:  jump of the solution at tₙ₊₁\nθ:  one-form evaluated across at tₙ\nΘ̅:  one-form evaluated across at tₙ₊₁\ng:  projection evaluated across  at tₙ\ng̅:  projection evaluated across at tₙ₊₁\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.NonlinearFunctionCacheDGVIPI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionCacheDGVIPI",
    "category": "type",
    "text": "Nonlinear function cache for Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nX: degrees of freedom\nQ: solution at quadrature nodes\nV: velocity at quadrature nodes\nP: one-form at quadrature nodes\nF: forces at quadrature nodes\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nq⁺: current solution of qₙ⁺\nq̅:  current solution of qₙ₊₁\nq̅⁻: current solution of qₙ₊₁⁻\nq̅⁺: current solution of qₙ₊₁⁺\nλ:  jump of the solution at tₙ\nλ̅:  jump of the solution at tₙ₊₁\nϕ:  solution evaluated across the jump at tₙ\nϕ̅:  solution evaluated across the jump at tₙ₊₁\nθ:  one-form evaluated across the jump at tₙ\nΘ̅:  one-form evaluated across the jump at tₙ₊₁\ng:  projection evaluated across the jump at tₙ\ng̅:  projection evaluated across the jump at tₙ₊₁\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.NonlinearFunctionCacheSIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionCacheSIPRK",
    "category": "type",
    "text": "Structure for holding the internal stages Q, the values of the drift vector and the diffusion matrix evaluated at the internal stages VQ=v(Q), BQ=B(Q), and the increments Y = Δta_driftv(Q) + a_diffB(Q)ΔW\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.NonlinearFunctionCacheSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionCacheSIRK",
    "category": "type",
    "text": "Structure for holding the internal stages Q, the values of the drift vector and the diffusion matrix evaluated at the internal stages VQ=v(Q), BQ=B(Q), and the increments Y = Δta_driftv(Q) + a_diffB(Q)ΔW\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.NonlinearFunctionCacheSISPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionCacheSISPRK",
    "category": "type",
    "text": "Structure for holding the internal stages Q, the values of the drift vector and the diffusion matrix evaluated at the internal stages VQ=v(Q), BQ=B(Q), and the increments Y = Δta_driftv(Q) + a_diffB(Q)ΔW\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersCGVI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersCGVI",
    "category": "type",
    "text": "ParametersCGVI: Parameters for right-hand side function of continuous Galerkin variational Integrator.\n\nParameters\n\nΘ: function of the noncanonical one-form (∂L/∂v)\nf: function of the force (∂L/∂q)\nΔt: time step\nb: weights of the quadrature rule\nc: nodes of the quadrature rule\nx: nodes of the basis\nm: mass matrix\na: derivative matrix\nr₀: reconstruction coefficients at the beginning of the interval\nr₁: reconstruction coefficients at the end of the interval\nt: current time\nq: current solution of q\np: current solution of p\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersDGVI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersDGVI",
    "category": "type",
    "text": "ParametersDGVI: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nDT: data type\nTT: parameter type\nD: dimension of the system\nS: number of basis nodes\nR: number of quadrature nodes\n\nFields\n\nΘ:  function of the noncanonical one-form (∂L/∂v)\nf:  function of the force (∂L/∂q)\ng:  function of the projection ∇ϑ(q)⋅v\nΔt: time step\nb:  quadrature weights\nc:  quadrature nodes\nm:  mass matrix\na:  derivative matrix\nr⁻: reconstruction coefficients, jump lhs value\nr⁺: reconstruction coefficients, jump rhs value\nt:  current time\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nq⁺: current solution of qₙ⁺\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersDGVIEXP",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersDGVIEXP",
    "category": "type",
    "text": "ParametersDGVIEXP: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nDT: data type\nTT: parameter type\nD: dimension of the system\nS: number of basis nodes\nR: number of quadrature nodes\n\nFields\n\nΘ:  function of the noncanonical one-form (∂L/∂v)\nf:  function of the force (∂L/∂q)\ng:  function of the projection ∇ϑ(q)⋅v\nΔt: time step\nb:  quadrature weights\nc:  quadrature nodes\nm:  mass matrix\na:  derivative matrix\nr⁻: reconstruction coefficients, jump lhs value\nr⁺: reconstruction coefficients, jump rhs value\nt:  current time\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nq⁺: current solution of qₙ⁺\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersDGVIP0",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersDGVIP0",
    "category": "type",
    "text": "ParametersDGVIP0: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nDT: data type\nTT: parameter type\nD: dimension of the system\nS: number of basis nodes\nR: number of quadrature nodes\n\nFields\n\nΘ:  function of the noncanonical one-form (∂L/∂v)\nf:  function of the force (∂L/∂q)\ng:  function of the projection ∇ϑ(q)⋅v\nΔt: time step\nb:  quadrature weights\nc:  quadrature nodes\nm:  mass matrix\na:  derivative matrix\nr⁻: reconstruction coefficients, jump lhs value\nr⁺: reconstruction coefficients, jump rhs value\nt:  current time\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nθ:  one-form ϑ evaluated on qₙ\nθ⁻: one-form ϑ evaluated on qₙ⁻\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersDGVIP1",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersDGVIP1",
    "category": "type",
    "text": "ParametersDGVIP1: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator.\n\nParameters\n\nDT: data type\nTT: parameter type\nD: dimension of the system\nS: number of basis nodes\nR: number of quadrature nodes\n\nFields\n\nΘ:  function of the noncanonical one-form (∂L/∂v)\nf:  function of the force (∂L/∂q)\ng:  function of the projection ∇ϑ(q)⋅v\nΔt: time step\nb:  quadrature weights\nc:  quadrature nodes\nm:  mass matrix\na:  derivative matrix\nr⁻: reconstruction coefficients, jump lhs value\nr⁺: reconstruction coefficients, jump rhs value\nt:  current time\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nq⁺: current solution of qₙ⁺\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersDGVIPI",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersDGVIPI",
    "category": "type",
    "text": "ParametersDGVIPI: Parameters for right-hand side function of Discontinuous Galerkin Variational Integrator with Path Integral approximation of the jump.\n\nParameters\n\nDT: data type\nTT: parameter type\nD:  dimension of the system\nS:  number of basis nodes\nQR: number of quadrature nodes\nFR: number of quadrature nodes for the discontinuity\n\nFields\n\nΘ:  function of the noncanonical one-form (∂L/∂v)\nf:  function of the force (∂L/∂q)\ng:  function of the projection ∇ϑ(q)⋅v\nΔt: time step\nb:  quadrature weights\nc:  quadrature nodes\nm:  mass matrix\na:  derivative matrix\nr⁻: reconstruction coefficients, jump lhs value\nr⁺: reconstruction coefficients, jump rhs value\nβ:  weights of the quadrature rule for the discontinuity\nγ:  nodes of the quadrature rule for the discontinuity\nμ⁻: mass vector for the lhs jump value\nμ⁺: mass vector for the rhs jump value\nα⁻: derivative vector for the discontinuity lhs value\nα⁺: derivative vector for the discontinuity rhs value\nρ⁻: reconstruction coefficients for central jump value, lhs value\nρ⁺: reconstruction coefficients for central jump value, rhs value\nt:  current time\nq:  current solution of qₙ\nq⁻: current solution of qₙ⁻\nq⁺: current solution of qₙ⁺\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersDIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersDIRK",
    "category": "type",
    "text": "Parameters for right-hand side function of diagonally implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersFIRK",
    "category": "type",
    "text": "Parameters for right-hand side function of fully implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersFLRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersFLRK",
    "category": "type",
    "text": "Parameters for right-hand side function of formal Lagrangian Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersHSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersHSPARK",
    "category": "type",
    "text": "Parameters for right-hand side function of Hamiltonian Specialised Partitioned Additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersIPRK",
    "category": "type",
    "text": "Parameters for right-hand side function of implicit partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersPARK",
    "category": "type",
    "text": "Parameters for right-hand side function of partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersPGLRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersPGLRK",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersSIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersSIPRK",
    "category": "type",
    "text": "Parameters for right-hand side function of implicit Runge-Kutta methods.   A - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersSIRK",
    "category": "type",
    "text": "Parameters for right-hand side function of implicit Runge-Kutta methods.    A - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersSISPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersSISPRK",
    "category": "type",
    "text": "Parameters for right-hand side function of implicit Runge-Kutta methods.   A - if positive, the upper bound of the Wiener process increments; if A=0.0, no truncation\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersSPARK",
    "category": "type",
    "text": "Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPARK",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRK",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKdegenerate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKdegenerate",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpInternal",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpInternal",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpLegendre",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpLegendre",
    "category": "type",
    "text": "Parameters for right-hand side function of variational special partitioned additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpMidpoint",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpMidpoint",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpSecondary",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpSecondary",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpStandard",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpStandard",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpSymmetric",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpSymmetric",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVPRKpVariational",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVPRKpVariational",
    "category": "type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersVSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersVSPARK",
    "category": "type",
    "text": "Parameters for right-hand side function of Variational Specialised Partitioned Additive Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.ParametersWIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.ParametersWIRK",
    "category": "type",
    "text": "Parameters for right-hand side function of weak implicit Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#Base.show-Tuple{IO,CoefficientsARK}",
    "page": "Integrators",
    "title": "Base.show",
    "category": "method",
    "text": "Print additive Runge-Kutta coefficients.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#Base.show-Tuple{IO,CoefficientsMRK}",
    "page": "Integrators",
    "title": "Base.show",
    "category": "method",
    "text": "Print multiplier Runge-Kutta coefficients.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#Base.show-Tuple{IO,CoefficientsPGLRK}",
    "page": "Integrators",
    "title": "Base.show",
    "category": "method",
    "text": "Print Runge-Kutta coefficients.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#Base.show-Tuple{IO,CoefficientsPRK}",
    "page": "Integrators",
    "title": "Base.show",
    "category": "method",
    "text": "Print projective Runge-Kutta coefficients.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#Base.show-Tuple{IO,CoefficientsRK}",
    "page": "Integrators",
    "title": "Base.show",
    "category": "method",
    "text": "Print Runge-Kutta coefficients.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.aitken_neville-Union{Tuple{TT}, Tuple{DT}, Tuple{Array{TT,1},Array{DT,2},TT,Array{DT,1}}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.aitken_neville",
    "category": "method",
    "text": "Compute p(x) where p is the unique polynomial of degree length(xi), such that p(x[i]) = y[i]) for all i.\n\nti: interpolation nodes\nxi: interpolation values\nt:  evaluation point\nx:  evaluation value\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.computeStageP!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorEPRK{DT,TT,VT,FT} where FT where VT,Int64,Int64,Int64,Any}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.computeStageP!",
    "category": "method",
    "text": "Compute P stages of explicit partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.computeStageQ!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorEPRK{DT,TT,VT,FT} where FT where VT,Int64,Int64,Int64,Any}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.computeStageQ!",
    "category": "method",
    "text": "Compute Q stages of explicit partitioned Runge-Kutta methods.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,2},1},Array{Array{ST,2},1},Array{Array{ST,2},1},Array{Array{ST,1},1},Array{Array{ST,1},1},ParametersSISPRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages!",
    "category": "method",
    "text": "Unpacks the data stored in x = (Y[1,1], Y[2,1], ... Y[D,1], Y[1,2], ..., Z[1,1], Z[2,1], ... Z[D,1], Z[1,2], ...) into the matrix Y, Z, calculates the internal stages Q, P, the values of the RHS of the SDE ( vi(Q,P), fi(Q,P), Bi(Q,P) and Gi(Q,P) ), and assigns them to VQPi, FQPi, BQPi and GQPi. Unlike for FIRK, here Y = Δt adrift v(Q,P) + adiff B(Q,P) ΔW, Z = Δt ̃a1drift f1(Q,P) + Δt ̃a2drift f2(Q,P) + ̃a1diff G1(Q,P) ΔW + ̃a2diff G2(Q,P) ΔW.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,2},1},Array{Array{ST,2},1},Array{Array{ST,1},1},Array{Array{ST,1},1},ParametersSIPRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages!",
    "category": "method",
    "text": "Unpacks the data stored in x = (Y[1,1], Y[2,1], ... Y[D,1], Y[1,2], ..., Z[1,1], Z[2,1], ... Z[D,1], Z[1,2], ...) into the matrix Y, Z, calculates the internal stages Q, P, the values of the RHS of the SDE ( v(Q,P), f(Q,P), B(Q,P) and G(Q,P) ), and assigns them to VQP, FQP, BQP and GQP. Unlike for FIRK, here Y = Δt adrift v(Q,P) + adiff B(Q,P) ΔW, Z = Δt ̃adrift v(Q,P) + ̃adiff B(Q,P) ΔW.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{Array{ST,1},1},Array{Array{ST,1},1},Array{Array{ST,2},1},Array{Array{ST,1},1},ParametersSIRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages!",
    "category": "method",
    "text": "Unpacks the data stored in x = (Y[1,1], Y[2,1], ... Y[D,1], Y[1,2], ...) into the matrix Y, calculates the internal stages Q, the values of the RHS of the SDE ( v(Q) and B(Q) ), and assigns them to VQ and BQ. Unlike for FIRK, here Y = Δt a v(Q) + ̃a B(Q) ΔW\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages!-Union{Tuple{S}, Tuple{M}, Tuple{D}, Tuple{ET}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{Array{ST,1},Array{Array{ST,1},1},Array{Array{ST,2},1},Array{Array{ST,1},1},Array{Array{ST,2},1},Array{Array{ST,1},1},Array{Array{ST,2},1},ParametersWIRK{DT,TT,ET,D,M,S}}} where S where M where D where ET where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages!",
    "category": "method",
    "text": "Unpacks the data stored in x = (Y0[1,1], Y0[2,1], ... Y0[D,1], ... Y0[D,S], Y1[1,1,1], Y1[2,1,1], ... Y1[D,1,1], Y1[1,2,1], Y1[2,2,1], ... Y1[D,2,1], ... Y1[D,M,S]  ) into the matrices Y0 and Y1, calculates the internal stages Q0 and Q1, the values of the RHS of the SDE ( v(Q0) and B(Q1) ), and assigns them to VQ and BQ. Unlike for FIRK, here Y = Δt a v(Q) + ̃a B(Q) ΔW\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_p!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVIPI{ST,D,S,R,FR} where FR,ParametersDGVIPI{DT,TT,D,S,R,FR,ΘT,FT,GT} where GT where FT where ΘT where FR}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_p!",
    "category": "method",
    "text": "Compute one-form and forces at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_p!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIEXP{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_p!",
    "category": "method",
    "text": "Compute one-form and forces at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_p!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIP0{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_p!",
    "category": "method",
    "text": "Compute one-form and forces at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_p!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIP1{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_p!",
    "category": "method",
    "text": "Compute one-form and forces at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_p!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVI{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_p!",
    "category": "method",
    "text": "Compute one-form and forces at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_q!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVIPI{ST,D,S,R,FR} where FR,ParametersDGVIPI{DT,TT,D,S,R,FR,ΘT,FT,GT} where GT where FT where ΘT where FR}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_q!",
    "category": "method",
    "text": "Compute solution at quadrature nodes and across jump.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_q!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIEXP{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_q!",
    "category": "method",
    "text": "Compute solution at quadrature nodes and across jump.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_q!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIP0{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_q!",
    "category": "method",
    "text": "Compute solution at quadrature nodes and across jump.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_q!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIP1{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_q!",
    "category": "method",
    "text": "Compute solution at quadrature nodes and across jump.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_q!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVI{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_q!",
    "category": "method",
    "text": "Compute solution at quadrature nodes and across jump.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_v!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVIPI{ST,D,S,R,FR} where FR,ParametersDGVIPI{DT,TT,D,S,R,FR,ΘT,FT,GT} where GT where FT where ΘT where FR}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_v!",
    "category": "method",
    "text": "Compute velocities at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_v!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIEXP{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_v!",
    "category": "method",
    "text": "Compute velocities at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_v!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIP0{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_v!",
    "category": "method",
    "text": "Compute velocities at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_v!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVIP1{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_v!",
    "category": "method",
    "text": "Compute velocities at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.compute_stages_v!-Union{Tuple{R}, Tuple{S}, Tuple{D}, Tuple{TT}, Tuple{DT}, Tuple{ST}, Tuple{NonlinearFunctionCacheDGVI{ST,D,S,R},ParametersDGVI{DT,TT,D,S,R,ΘT,FT,GT} where GT where FT where ΘT}} where R where S where D where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.compute_stages_v!",
    "category": "method",
    "text": "Compute velocities at quadrature nodes.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.create_internal_stage_vector-NTuple{4,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.create_internal_stage_vector",
    "category": "method",
    "text": "Create a vector of S solution matrices of type DT to store the solution of S internal stages for a problem with DxM dimensions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.create_internal_stage_vector-Tuple{Any,Any,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.create_internal_stage_vector",
    "category": "method",
    "text": "Create a vector of S solution vectors of type DT to store the solution of S internal stages for a problem with D dimensions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.create_internal_stage_vector_with_zero-Tuple{Any,Any,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.create_internal_stage_vector_with_zero",
    "category": "method",
    "text": "Create a vector of S+1 solution vectors of type DT to store the solution of S internal stages and the solution of the previous timestep for a problem with D     dimensions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.create_nonlinear_solver-Tuple{Any,Any,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.create_nonlinear_solver",
    "category": "method",
    "text": "Create nonlinear solver object for a system of N equations with data type DT. The function f(x)=0 to be solved for is determined by a julia function function_stages!(x, b, params), where x is the current solution and b is the output vector, s.th. b = f(x). params are a set of parameters depending on the equation and integrator that is used. The solver type is obtained from the config dictionary (:nls_solver).\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.create_solution_vector-NTuple{4,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.create_solution_vector",
    "category": "method",
    "text": "Create a solution vector of type TwicePrecision{DT} for a problem with D dimensions, NS\' sample paths, andNI` independent initial conditions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.create_solution_vector-Tuple{Any,Any,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.create_solution_vector",
    "category": "method",
    "text": "Create a solution vector of type TwicePrecision{DT} for a problem with D dimensions and M independent initial conditions.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.euler_extrapolation-Union{Tuple{TT}, Tuple{DT}, Tuple{Function,TT,TT,Array{DT,1},Array{DT,1},Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.euler_extrapolation",
    "category": "method",
    "text": "Euler extrapolation method with arbitrary order p.\n\nv:  function to compute vector field\nt₀: initial time\nt₁: final   time\nx₀: initial value\nx₁: final   value\ns:  number of interpolations (order p=s+1)\n\nTODO This is probably broken!\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initial_guess!-Tuple{IntegratorDIRK,Int64}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initial_guess!",
    "category": "method",
    "text": "Compute initial guess for internal stages.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initial_guess!-Union{Tuple{IntegratorSIPRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersSIPRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(PSDE{DT,TT,vType,fType,BType,GType,N} where N where GType<:Function where BType<:Function where fType<:Function where vType<:Function))}, Tuple{TT}, Tuple{DT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initial_guess!",
    "category": "method",
    "text": "This function computes initial guesses for Y, Z and assigns them to int.solver.x The prediction is calculated using an explicit integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initial_guess!-Union{Tuple{IntegratorSIRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersSIRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(SDE{DT,TT,vType,BType,N} where N where BType<:Function where vType<:Function))}, Tuple{TT}, Tuple{DT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initial_guess!",
    "category": "method",
    "text": "This function computes initial guesses for Y and assigns them to int.solver.x The prediction is calculated using an explicit integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initial_guess!-Union{Tuple{IntegratorSISPRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersSISPRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(SPSDE{DT,TT,vType,f1Type,f2Type,BType,G1Type,G2Type,N} where N where G2Type<:Function where G1Type<:Function where BType<:Function where f2Type<:Function where f1Type<:Function where vType<:Function))}, Tuple{TT}, Tuple{DT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initial_guess!",
    "category": "method",
    "text": "This function computes initial guesses for Y, Z and assigns them to int.solver.x The prediction is calculated using an explicit integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.initial_guess!-Union{Tuple{IntegratorWIRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersWIRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(SDE{DT,TT,vType,BType,N} where N where BType<:Function where vType<:Function))}, Tuple{TT}, Tuple{DT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initial_guess!",
    "category": "method",
    "text": "This function computes initial guesses for Y and assigns them to int.solver.x The prediction is calculated using an explicit integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{NW}, Tuple{NQ}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorSERK{DT,TT,FT},SolutionSDE{DT,TT,NQ,NW},Int64,Int64,Int64}} where NW where NQ where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate SDE with explicit Runge-Kutta integrator.   Calculating the n-th time step of the explicit integrator for the   sample path r and the initial condition m\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{NW}, Tuple{NQ}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorWERK{DT,TT,FT},SolutionSDE{DT,TT,NQ,NW},Int64,Int64,Int64}} where NW where NQ where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate SDE with explicit Runge-Kutta integrator.  Calculating the n-th time step of the explicit integrator for the sample path r and the initial condition m\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorSIPRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersSIPRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(PSDE{DT,TT,vType,fType,BType,GType,N} where N where GType<:Function where BType<:Function where fType<:Function where vType<:Function)),SolutionPSDE{DT,TT,NQ,NW},Int64,Int64,Int64}} where NW where NQ where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate PSDE with a stochastic implicit partitioned Runge-Kutta integrator. Integrating the k-th sample path for the m-th initial condition\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorSIRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersSIRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(SDE{DT,TT,vType,BType,N} where N where BType<:Function where vType<:Function)),SolutionSDE{DT,TT,NQ,NW},Int64,Int64,Int64}} where NW where NQ where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate SDE with a stochastic implicit Runge-Kutta integrator.   Integrating the k-th sample path for the m-th initial condition\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorSISPRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersSISPRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(SPSDE{DT,TT,vType,f1Type,f2Type,BType,G1Type,G2Type,N} where N where G2Type<:Function where G1Type<:Function where BType<:Function where f2Type<:Function where f1Type<:Function where vType<:Function)),SolutionPSDE{DT,TT,NQ,NW},Int64,Int64,Int64}} where NW where NQ where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate PSDE with a stochastic implicit partitioned Runge-Kutta integrator.  Integrating the k-th sample path for the m-th initial condition\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorWIRK{DT,TT,PT,ST,N} where N where ST<:NonlinearSolver{DT} where PT<:(ParametersWIRK{DT,TT,ET,D,M,S} where S where M where D where ET<:(SDE{DT,TT,vType,BType,N} where N where BType<:Function where vType<:Function)),SolutionSDE{DT,TT,NQ,NW},Int64,Int64,Int64}} where NW where NQ where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate SDE with a stochastic implicit Runge-Kutta integrator.   Integrating the k-th sample path for the m-th initial condition\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorSplitting{DT,TT,FT,ST,FT1,CT,N} where N where CT where FT1 where ST<:AbstractTableauSplitting,SolutionODE{DT,TT,N},Int64,Int64}} where N where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with splitting integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorDIRK{DT,TT,PT,ST,IT,N} where N where IT<:(InitialGuessODE{DT,TT,VT,IT} where IT<:Interpolator where VT) where ST where PT<:(ParametersDIRK{DT,TT,ET,D,S} where S where D where ET<:(ODE{DT,TT,vType,N} where N where vType<:Function)),SolutionODE{DT,TT,N},Int64,Int64}} where N where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with diagonally implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorFIRK{DT,TT,PT,ST,IT,N} where N where IT<:(InitialGuessODE{DT,TT,VT,IT} where IT<:Interpolator where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersFIRK{DT,TT,ET,D,S} where S where D where ET<:(ODE{DT,TT,vType,N} where N where vType<:Function)),SolutionODE{DT,TT,N},Int64,Int64}} where N where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with fully implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{VT}, Tuple{FT}, Tuple{ΘT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpLegendre{DT,TT,ΘT,FT,VT,VT1,SPT,ST,IT} where IT where ST where SPT where VT1,SolutionPDAE{DT,TT,N},Int64,Int64}} where N where VT where FT where ΘT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate DAE with variational special partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{VT}, Tuple{ϕT}, Tuple{GT}, Tuple{UT}, Tuple{PT}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorSPARK{DT,TT,FT,PT,UT,GT,ϕT,VT,SPT,ST,IT} where IT where ST where SPT,SolutionPDAE{DT,TT,N},Int64,Int64}} where N where VT where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate DAE with variational special partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{VT}, Tuple{ϕT}, Tuple{GT}, Tuple{UT}, Tuple{PT}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPARK{DT,TT,FT,PT,UT,GT,ϕT,VT,SPT,ST,IT} where IT where ST where SPT,SolutionPDAE{DT,TT,N},Int64,Int64}} where N where VT where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate DAE with variational partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{ψT}, Tuple{ϕT}, Tuple{FT}, Tuple{VT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorHSPARK{DT,TT,VT,FT,ϕT,ψT,SPT,ST,IT} where IT where ST where SPT,SolutionPDAE{DT,TT,N},Int64,Int64}} where N where ψT where ϕT where FT where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate DAE with variational special partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{N}, Tuple{ϕT}, Tuple{GT}, Tuple{UT}, Tuple{PT}, Tuple{FT}, Tuple{TT}, Tuple{DT}, Tuple{IntegratorPARK{DT,TT,FT,PT,UT,GT,ϕT,ST} where ST,SolutionPDAE{DT,TT,N},Int64,Int64}} where N where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate DAE with partitioned additive Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorCGVI{DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT,D,S,R} where R where S where D where BT<:Basis where IT where ST where FPT where VT where GT where FT where ΘT,Union{SolutionPDAE{DT,TT,N} where N, SolutionPODE{DT,TT,N} where N},Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorDGVIEXP{DT,TT,D,S,R,ΘT,FT,GT,VT,FPT,ST,IT,BT} where BT<:Basis where IT where ST where FPT where VT where GT where FT where ΘT where R where S where D,Union{SolutionPDAE{DT,TT,N} where N, SolutionPODE{DT,TT,N} where N},Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorDGVIP0{DT,TT,D,S,R,ΘT,FT,GT,VT,FPT,ST,IT,BT} where BT<:Basis where IT where ST where FPT where VT where GT where FT where ΘT where R where S where D,Union{SolutionPDAE{DT,TT,N} where N, SolutionPODE{DT,TT,N} where N},Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorDGVIP1{DT,TT,D,S,R,ΘT,FT,GT,VT,FPT,ST,IT,BT} where BT<:Basis where IT where ST where FPT where VT where GT where FT where ΘT where R where S where D,Union{SolutionPDAE{DT,TT,N} where N, SolutionPODE{DT,TT,N} where N},Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorDGVIPI{DT,TT,D,S,R,ΘT,FT,GT,VT,FPT,ST,IT,BT,JT} where JT<:Discontinuity where BT<:Basis where IT where ST where FPT where VT where GT where FT where ΘT where R where S where D,Union{SolutionPDAE{DT,TT,N} where N, SolutionPODE{DT,TT,N} where N},Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorDGVI{DT,TT,D,S,R,ΘT,FT,GT,VT,FPT,ST,IT,BT} where BT<:Basis where IT where ST where FPT where VT where GT where FT where ΘT where R where S where D,Union{SolutionPDAE{DT,TT,N} where N, SolutionPODE{DT,TT,N} where N},Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorEPRK{DT,TT,VT,FT} where FT where VT,SolutionPODE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorERK{DT,TT,FT} where FT,SolutionODE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with explicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorFLRK{DT,TT,AT,FT,GT,VT,ΩT,dHT,SPT,ST,IT,N} where N where IT<:(InitialGuessODE{DT,TT,VT,IT} where IT<:Interpolator) where ST where SPT where dHT where ΩT where VT where GT where FT where AT,SolutionPODE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with fully implicit Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorIPRK{DT,TT,PT,ST,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersIPRK{DT,TT,ET,D,S} where S where D where ET<:(PODE{DT,TT,vType,fType,N} where N where fType<:Function where vType<:Function)),SolutionPODE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with implicit partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorPGLRK{DT,TT,PT,ST,IT,N} where N where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersPGLRK{DT,TT,D,S,ET} where ET where S where D),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate PODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKdegenerate{DT,TT,SPT,PPT,SST,STP,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where STP<:NonlinearSolver{DT} where SST<:NonlinearSolver{DT} where PPT<:(ParametersVPRKdegenerate{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)) where SPT<:(ParametersVPRK{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpInternal{DT,TT,PT,ST,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersVPRKpInternal{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpMidpoint{DT,TT,PT,ST,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersVPRKpMidpoint{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpSecondary{DT,TT,PT,ST,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersVPRKpSecondary{DT,TT,ET,D,S} where S where D where ET<:(VODE{DT,TT,αType,fType,gType,vType,ωType,dHType,N} where N where dHType<:Function where ωType<:Function where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpStandard{DT,TT,SPT,PPT,SST,STP,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where STP<:NonlinearSolver{DT} where SST<:NonlinearSolver{DT} where PPT<:(ParametersVPRKpStandard{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)) where SPT<:(ParametersVPRK{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpSymmetric{DT,TT,PT,ST,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersVPRKpSymmetric{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRKpVariational{DT,TT,SPT,PPT,SST,STP,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where STP<:NonlinearSolver{DT} where SST<:NonlinearSolver{DT} where PPT<:(ParametersVPRKpVariational{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)) where SPT<:(ParametersVPRK{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{TT}, Tuple{DT}, Tuple{IntegratorVPRK{DT,TT,PT,ST,IT} where IT<:(InitialGuessPODE{DT,TT,VT,FT,IT} where IT<:Interpolator where FT where VT) where ST<:NonlinearSolver{DT} where PT<:(ParametersVPRK{DT,TT,ET,D,S} where S where D where ET<:(IODE{DT,TT,αType,fType,gType,vType,N} where N where vType<:Function where gType<:Function where fType<:Function where αType<:Function)),SolutionPDAE{DT,TT,N} where N,Int64,Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.midpoint_extrapolation-Union{Tuple{TT}, Tuple{DT}, Tuple{Function,Function,TT,TT,Array{DT,1},Array{DT,1},Array{DT,1},Array{DT,1},Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.midpoint_extrapolation",
    "category": "method",
    "text": "Midpoint extrapolation method with arbitrary order p.\n\nv:  function to compute vector field\nf:  function to compute force  field\nt₀: initial time\nt₁: final   time\nq₀: initial positions\np₀: initial momenta\nq₁: final   positions\np₁: final   momenta\ns:  number of interpolations (order p=2s+2)\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.midpoint_extrapolation-Union{Tuple{TT}, Tuple{DT}, Tuple{Function,TT,TT,Array{DT,1},Array{DT,1},Int64}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.midpoint_extrapolation",
    "category": "method",
    "text": "Midpoint extrapolation method with arbitrary order p.\n\nv:  function to compute vector field\nt₀: initial time\nt₁: final   time\nx₀: initial value\nx₁: final   value\ns:  number of interpolations (order p=2s+2)\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#GeometricIntegrators.Integrators.readTableauRKHeaderFromFile-Tuple{Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.readTableauRKHeaderFromFile",
    "category": "method",
    "text": "Reads and parses Tableau metadata from file.\n\n\n\n\n\n"
},

{
    "location": "modules/integrators/#Integrators-1",
    "page": "Integrators",
    "title": "Integrators",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Integrators]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/interpolation/#",
    "page": "Interpolation",
    "title": "Interpolation",
    "category": "page",
    "text": ""
},

{
    "location": "modules/interpolation/#GeometricIntegrators.Interpolation.HermiteInterpolation",
    "page": "Interpolation",
    "title": "GeometricIntegrators.Interpolation.HermiteInterpolation",
    "category": "type",
    "text": "Hermite\'s Interpolating Polynomials\n\nHere, we implement a two point Hermite interpolation function which passes through the function and its first derivative for the interval 01. The polynomial is determined by four constraint equations, matching the function and its derivative at the points 0 and 1.\n\nStart by defining the 3rd degree polynomial and its derivative by\n\nbeginalign*\ng(x) = a_0 + a_1 x + a_2 x^2 + a_3 x^3  \ng(x) = a_1 + 2 a_2 x + 3 a_3 x^2 \nendalign*\n\nand apply the constraints\n\nbeginalign*\ng(0) = f_0   Rightarrow  a_0 = f_0  \ng(1) = f_1   Rightarrow  a_0 + a_1 + a_2 + a_3 = f_1  \ng(0) = f_0   Rightarrow  a_1 = f_0  \ng(1) = f_1   Rightarrow  a_1 + 2 a_2 + 3 a_3 = f_1  \nendalign*\n\nSolving for a_0 a_1 a_2 a_3 leads to\n\nbeginalign*\na_0 = f_0  \na_1 = f_0  \na_2 = - 3 f_0 + 3 f_1 - 2 f_0 - f_1  \na_3 = 2 f_0 - 2 f_1 + f_0 + f_1 \nendalign*\n\nso that the polynomial g(x) reads\n\ng(x) = f_0 + f_0 x + (- 3 f_0 + 3 f_1 - 2 f_0 - f_1) x^2 + (2 f_0 - 2 f_1 + f_0 + f_1) x^3 \n\nThe function and derivative values can be factored out, so that g(x) can be rewritten as\n\ng(x) = f_0 (1 - 3 x^2 + 2 x^3) + f_1 (3 x^2 - 2 x^3) + f_0 (x - 2 x^2 + x^3) + f_1 (- x^2 + x^3) \n\nor in generic form as\n\ng(x) = f_0 a_0(x) + f_1 a_1(x) + f_0 b_0(x) + f_1 b_1(x) \n\nwith basis functions\n\nbeginalign*\na_0 (x) = 1 - 3 x^2 + 2 x^3  \nb_0 (x) = x - 2 x^2 + x^3  \na_1 (x) = 3 x^2 - 2 x^3  \nb_1 (x) = - x^2 + x^3 \nendalign*\n\nThe derivative g(x) accordingly reads\n\ng(x) = f_0 a_0(x) + f_1 a_1(x) + f_0 b_0(x) + f_1 b_1(x) \n\nwith\n\nbeginalign*\na_0 (x) = - 6 x + 6 x^2  \nb_0 (x) = 1 - 4 x + 3 x^2  \na_1 (x) = 6 x - 6 x^2  \nb_1 (x) = - 2 x + 3 x^2 \nendalign*\n\nThe basis functions a_0and a_1 are associated with the function values at x_0 and x_1, respectively, while the basis functions b_0 and b_1 are associated with the derivative values at x_0 and x_1. The basis functions satisfy the following relations,\n\nbeginalign*\na_i (x_j) = delta_ij  \nb_i (x_j) = 0  \na_i (x_j) = 0  \nb_i (x_j) = delta_ij  \nij = 0 1 \nendalign*\n\nwhere delta_ij denotes the Kronecker-delta, so that\n\nbeginalign*\ng(0) = f_0  \ng(1) = f_1  \ng(0) = f_0  \ng(1) = f_1 \nendalign*\n\n\n\n\n\n"
},

{
    "location": "modules/interpolation/#Interpolation-1",
    "page": "Interpolation",
    "title": "Interpolation",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Interpolation]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solvers_linear/#",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solvers_linear/#Linear-Solvers-1",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solvers]\nPages   = [\"solvers/linear/linear_solvers.jl\",\n           \"solvers/linear/lu_solver_lapack.jl\"]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solvers_nonlinear/#",
    "page": "Nonlinear Solvers",
    "title": "Nonlinear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solvers_nonlinear/#Nonlinear-Solvers-1",
    "page": "Nonlinear Solvers",
    "title": "Nonlinear Solvers",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solvers]\nPages   = [\"solvers/nonlinear/nonlinear_solvers.jl\",\n           \"solvers/nonlinear/jacobian.jl\",\n           \"solvers/nonlinear/abstract_fixed_point_solver.jl\",\n           \"solvers/nonlinear/fixed_point_solver.jl\",\n           \"solvers/nonlinear/abstract_newton_solver.jl\",\n           \"solvers/nonlinear/newton_solver.jl\",\n           \"solvers/nonlinear/quasi_newton_solver.jl\"]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/quadratures/#",
    "page": "Quadrature Rules",
    "title": "Quadrature Rules",
    "category": "page",
    "text": ""
},

{
    "location": "modules/quadratures/#GeometricIntegrators.Quadratures.quadrature-Tuple{Quadrature,Function}",
    "page": "Quadrature Rules",
    "title": "GeometricIntegrators.Quadratures.quadrature",
    "category": "method",
    "text": "Integrate a function f(x) over the interval [0,1] using the quadrature quad.\n\n\n\n\n\n"
},

{
    "location": "modules/quadratures/#GeometricIntegrators.Quadratures.shift!-Tuple{Any,Any}",
    "page": "Quadrature Rules",
    "title": "GeometricIntegrators.Quadratures.shift!",
    "category": "method",
    "text": "Scale nodes and weights from the interval [-1,+1] to the interval [0,1].\n\n\n\n\n\n"
},

{
    "location": "modules/quadratures/#GeometricIntegrators.Quadratures.unshift!-Tuple{Any,Any}",
    "page": "Quadrature Rules",
    "title": "GeometricIntegrators.Quadratures.unshift!",
    "category": "method",
    "text": "Scale nodes and weights from the interval [0,1] to the interval [-1,+1].\n\n\n\n\n\n"
},

{
    "location": "modules/quadratures/#Quadratures-1",
    "page": "Quadrature Rules",
    "title": "Quadratures",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Quadratures]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solutions/#",
    "page": "Solutions",
    "title": "Solutions",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.PSolutionPDAE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.PSolutionPDAE",
    "category": "type",
    "text": "Parallel Solution of a partitioned differential algebraic equation.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.SSolutionPDAE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SSolutionPDAE",
    "category": "type",
    "text": "Serial Solution of a partitioned differential algebraic equation.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for partitioned DAE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for ODE and split ODE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for partitioned ODE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for SDE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for DAE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for variational ODE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for PSDE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Print error for solutions of equations not implemented, yet.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for implicit DAE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "type",
    "text": "Create solution for implicit ODE.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.SolutionDAE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionDAE",
    "category": "type",
    "text": "Solution of a differential algebraic equation.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.SolutionODE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionODE",
    "category": "type",
    "text": "SolutionODE: Solution of an ordinary differential equation\n\nContains all fields necessary to store the solution of an ODE.\n\nFields\n\nnd: dimension of the dynamical variable q\nnt: number of time steps to store\nni: number of initial conditions\nt:  time steps\nq:  solution q[nd, nt+1, ni] with q[:,0,:] the initial conditions\nntime: number of time steps to compute\nnsave: save every nsave\'th time step\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.SolutionPODE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionPODE",
    "category": "type",
    "text": "Solution of a partitioned ordinary differential equation.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.SolutionPSDE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionPSDE",
    "category": "type",
    "text": "SolutionPSDE: Solution of a partitioned stochastic differential equation\n\nContains all fields necessary to store the solution of a PSDE or SPSDE\n\nFields\n\nconv: type of the solution: \"strong\" or \"weak\"\nnd: dimension of the dynamical variable q\nnm: dimension of the Wiener process\nnt: number of time steps to store\nns: number of sample paths\nni: number of initial conditions\nt:  time steps\nq:  solution q[nd, nt+1, ns, ni] with q[:,0,:,:] the initial conditions\np:  solution p[nd, nt+1, ns, ni] with p[:,0,:,:] the initial conditions\nW:  Wiener process driving the stochastic processes q and p\nK:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions),\n  A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; if K=0 no truncation\nntime: number of time steps to compute\nnsave: save every nsave\'th time step\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.SolutionSDE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionSDE",
    "category": "type",
    "text": "SolutionSDE: Solution of a stochastic differential equation\n\nContains all fields necessary to store the solution of an SDE.\n\nFields\n\nconv: type of the solution: \"strong\" or \"weak\"\nnd: dimension of the dynamical variable q\nnm: dimension of the Wiener process\nnt: number of time steps to store\nns: number of sample paths\nni: number of initial conditions\nt:  time steps\nq:  solution q[nd, nt+1, ns, ni] with q[:,0,:,:] the initial conditions\nW:  Wiener process driving the stochastic process q\nK:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions),\n  A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; if K=0 no truncation\nntime: number of time steps to compute\nnsave: save every nsave\'th time step\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{SolutionPSDE{DT,TT,NQ,NW},HDF5File}, Tuple{SolutionPSDE{DT,TT,NQ,NW},HDF5File,Any}, Tuple{SolutionPSDE{DT,TT,NQ,NW},HDF5File,Any,Any}} where NW where NQ where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.   offset - start writing q at the position offset+2   offset2- start writing ΔW, ΔZ at the position offset2+1\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{SolutionSDE{DT,TT,NQ,NW},HDF5File}, Tuple{SolutionSDE{DT,TT,NQ,NW},HDF5File,Any}, Tuple{SolutionSDE{DT,TT,NQ,NW},HDF5File,Any,Any}} where NW where NQ where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.   offset - start writing q at the position offset+2   offset2- start writing ΔW, ΔZ at the position offset2+1\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionDAE{DT,TT,2},HDF5File}, Tuple{SolutionDAE{DT,TT,2},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionDAE{DT,TT,3},HDF5File}, Tuple{SolutionDAE{DT,TT,3},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionODE{DT,TT,2},HDF5File}, Tuple{SolutionODE{DT,TT,2},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionODE{DT,TT,3},HDF5File}, Tuple{SolutionODE{DT,TT,3},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPDAE{DT,TT,2},HDF5File}, Tuple{SolutionPDAE{DT,TT,2},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPDAE{DT,TT,3},HDF5File}, Tuple{SolutionPDAE{DT,TT,3},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPODE{DT,TT,2},HDF5File}, Tuple{SolutionPODE{DT,TT,2},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.CommonFunctions.write_to_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPODE{DT,TT,3},HDF5File}, Tuple{SolutionPODE{DT,TT,3},HDF5File,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "method",
    "text": "Append solution to HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{SolutionPSDE{DT,TT,NQ,NW},AbstractString}, Tuple{SolutionPSDE{DT,TT,NQ,NW},AbstractString,Int64}, Tuple{SolutionPSDE{DT,TT,NQ,NW},AbstractString,Int64,Int64}} where NW where NQ where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for SDE solution object.   It is implemented as one fucntion for all NQ and NW cases, rather than several   separate cases as was done for SolutionODE.   nt - the total number of time steps to store   ntime - the total number of timesteps to be computed\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{NW}, Tuple{NQ}, Tuple{TT}, Tuple{DT}, Tuple{SolutionSDE{DT,TT,NQ,NW},AbstractString}, Tuple{SolutionSDE{DT,TT,NQ,NW},AbstractString,Int64}, Tuple{SolutionSDE{DT,TT,NQ,NW},AbstractString,Int64,Int64}} where NW where NQ where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for SDE solution object.   It is implemented as one fucntion for all NQ and NW cases, rather than several   separate cases as was done for SolutionODE.   nt - the total number of time steps to store   ntime - the total number of timesteps to be computed\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionDAE{DT,TT,2},AbstractString}, Tuple{SolutionDAE{DT,TT,2},AbstractString,Int64}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for DAE solution object.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionDAE{DT,TT,3},AbstractString}, Tuple{SolutionDAE{DT,TT,3},AbstractString,Int64}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for DAE solution object.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionODE{DT,TT,2},AbstractString}, Tuple{SolutionODE{DT,TT,2},AbstractString,Int64}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for ODE solution object.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionODE{DT,TT,3},AbstractString}, Tuple{SolutionODE{DT,TT,3},AbstractString,Int64}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for ODE solution object.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPDAE{DT,TT,2},AbstractString}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for PDAE solution object with single initial condition.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPDAE{DT,TT,3},AbstractString}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for PDAE solution object with multiple initial conditions.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPODE{DT,TT,2},AbstractString}, Tuple{SolutionPODE{DT,TT,2},AbstractString,Int64}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for PODE solution object.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPODE{DT,TT,3},AbstractString}, Tuple{SolutionPODE{DT,TT,3},AbstractString,Int64}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "method",
    "text": "Creates HDF5 file and initialises datasets for PODE solution object.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.get_initial_conditions!-Union{Tuple{TT}, Tuple{DT}, Tuple{SolutionPSDE{DT,TT,NQ,NW} where NW where NQ,Union{Array{DT,1}, Array{TwicePrecision{DT},1}},Union{Array{DT,1}, Array{TwicePrecision{DT},1}},Any,Any}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.get_initial_conditions!",
    "category": "method",
    "text": "copies the m-th initial condition for the k-th sample path from sol.q to q\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.createHDF5",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.createHDF5",
    "category": "function",
    "text": "createHDF5: Creates or opens HDF5 file.   A version for StochasticSolution. It does not create attributes   and does not write the time array t, like the version above does. Instead these   are set in create_hdf5(), so that arrays larger than currently held in the solution   structure can be created in the file. In the future it would be better to rewrite   the function above, so that it is universal for all solution structures.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.createHDF5",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.createHDF5",
    "category": "function",
    "text": "createHDF5: Creates or opens HDF5 file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#GeometricIntegrators.Solutions.writeSolutionToHDF5-Tuple{Solution,AbstractString}",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.writeSolutionToHDF5",
    "category": "method",
    "text": "Creates HDF5 file, writes solution to file, and closes file.\n\n\n\n\n\n"
},

{
    "location": "modules/solutions/#Solutions-1",
    "page": "Solutions",
    "title": "Solutions",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solutions]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/tableaus/#",
    "page": "Tableaus",
    "title": "Tableaus",
    "category": "page",
    "text": ""
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauBurrageCL-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauBurrageCL",
    "category": "method",
    "text": "Tableau for the explicit 4-stage CL method due to K. Burrage and P. Burrage   Method cited in Eq. (56) in    K. Burrage, P. Burrage (1996) \"High strong order explicit Runge-Kutta methods for stochastic ordinary differential equations\".   According to the paper, the method has strong order 1.5 for one-dimensional Brownian motion. Reduces to the classical R-K method   of order 4 when noise is zero.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauBurrageE1-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauBurrageE1",
    "category": "method",
    "text": "Tableau for the explicit 4-stage E1 method due to K. Burrage and P. Burrage   Method cited in Eq. (4.2)-(4.3) in    K. Burrage, P. Burrage (2000) \"Order conditions for stochastic Runge-Kutta methods by B-series\".   According to the paper, the method has strong order 1.0 for one-dimensional Brownian motion.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauBurrageG5-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauBurrageG5",
    "category": "method",
    "text": "Tableau for the explicit 5-stage G5 method due to K. Burrage and P. Burrage   Method cited in Section 4 of    K. Burrage, P. Burrage (2000) \"Order conditions for stochastic Runge-Kutta methods by B-series\".   According to the paper, the method has strong order 1.5 for one-dimensional Brownian motion.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauBurrageR2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauBurrageR2",
    "category": "method",
    "text": "Tableau for the explicit 2-stage R2 method due to K. Burrage and P. Burrage   Method cited in Eq. (51) in    K. Burrage, P. Burrage (1996) \"High strong order explicit Runge-Kutta methods for stochastic ordinary differential equations\".   According to the paper, the method has strong order 1.0 for one-dimensional Brownian motion\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauERK4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauERK4",
    "category": "method",
    "text": "Tableau for explicit Runge-Kutta method of order four (1/6 rule)\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauERK438-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauERK438",
    "category": "method",
    "text": "Tableau for explicit Runge-Kutta method of order four (3/8 rule)\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauExplicitEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauExplicitEuler",
    "category": "method",
    "text": "Tableau for explicit Euler method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauExplicitMidpoint-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauExplicitMidpoint",
    "category": "method",
    "text": "Tableau for explicit midpoint method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauGLRKpSymmetric-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRKpSymmetric",
    "category": "method",
    "text": "Tableau for Gauss-Legendre method with s stages and symplectic projection.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauGLRKpSymplectic-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRKpSymplectic",
    "category": "method",
    "text": "Tableau for Gauss-Legendre method with s stages and symplectic projection.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauHeun-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauHeun",
    "category": "method",
    "text": "Tableau for Heun\'s method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauImplicitEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauImplicitEuler",
    "category": "method",
    "text": "Implicit Euler\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauImplicitMidpoint-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauImplicitMidpoint",
    "category": "method",
    "text": "Implicit Midpoint\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauKutta-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauKutta",
    "category": "method",
    "text": "Tableau for Kutta\'s method of order three\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIA2",
    "category": "method",
    "text": "Gauss-Lobatto-IIIA Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIA3",
    "category": "method",
    "text": "Gauss-Lobatto-IIIA Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIA4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIA4",
    "category": "method",
    "text": "Gauss-Lobatto-IIIA Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymmetric-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymmetric",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymplectic-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymplectic",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symplectic projection.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymmetric-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymmetric",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymplectic-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymplectic",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symplectic projection.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIB2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIB2",
    "category": "method",
    "text": "Gauss-Lobatto-IIIB Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIB3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIB3",
    "category": "method",
    "text": "Gauss-Lobatto-IIIB Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIB4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIB4",
    "category": "method",
    "text": "Gauss-Lobatto-IIIB Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIC2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIC2",
    "category": "method",
    "text": "Gauss-Lobatto-IIIC Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIC3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIC3",
    "category": "method",
    "text": "Gauss-Lobatto-IIIC Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIC4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIC4",
    "category": "method",
    "text": "Gauss-Lobatto-IIIC Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIID2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIID2",
    "category": "method",
    "text": "Gauss-Lobatto-IIID Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIID3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIID3",
    "category": "method",
    "text": "Gauss-Lobatto-IIID Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIID4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIID4",
    "category": "method",
    "text": "Gauss-Lobatto-IIID Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIE2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIE2",
    "category": "method",
    "text": "Gauss-Lobatto-IIIE Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIE3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIE3",
    "category": "method",
    "text": "Gauss-Lobatto-IIIE Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIE4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIE4",
    "category": "method",
    "text": "Gauss-Lobatto-IIIE Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIF2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIF2",
    "category": "method",
    "text": "Gauss-Lobatto-IIIF Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIF3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIF3",
    "category": "method",
    "text": "Gauss-Lobatto-IIIF Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIF4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIF4",
    "category": "method",
    "text": "Gauss-Lobatto-IIIF Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIG2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIG2",
    "category": "method",
    "text": "Gauss-Lobatto-IIIG Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIG3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIG3",
    "category": "method",
    "text": "Gauss-Lobatto-IIIG Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobIIIG4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIG4",
    "category": "method",
    "text": "Gauss-Lobatto-IIIG Runge-Kutta, s=4\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobattoIIIAIIIB2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobattoIIIAIIIB2",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIAIIIB method with s=2 stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauLobattoIIIBIIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobattoIIIBIIIA2",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIBIIIA method with s=2 stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauModifiedStochasticStormerVerlet",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauModifiedStochasticStormerVerlet",
    "category": "function",
    "text": "Tableau for the 2-stage modified stochastic LobattoIIIA-IIIB method   Tableau for the 2-stage modified stochastic LobattoIIIA-IIIB method   Satisfies the conditions for Lagrange-d\'Alembert integrators   and the conditions for convergence of order 1.0 for one Wiener process\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauPlaten-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauPlaten",
    "category": "method",
    "text": "Tableau for the explicit Platen method   Platen\'s method cited in Eq. (52) in    K. Burrage, P. Burrage (1996) \"High strong order explicit Runge-Kutta methods for stochastic ordinary differential equations\".   According to the paper, the method has strong order 1.0 for one-dimensional Brownian motion.   Appears to have a rather poor long-time performance.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauRadIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauRadIIA2",
    "category": "method",
    "text": "Gauss-Radau-IIA Runge-Kutta, s=2\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauRadIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauRadIIA3",
    "category": "method",
    "text": "Gauss-Radau-IIA Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauRosslerRS1-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauRosslerRS1",
    "category": "method",
    "text": "Tableau for the explicit 4-stage RS1 method due to Andreas Rossler   Method cited in Table 5.2 in   Andreas Rossler, \"Second order Runge-Kutta methods for Stratonovich stochastic differential equations\",   BIT Numerical Mathematics (2007) 47   According to the paper, the method has weak order 2.0.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauRosslerRS2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauRosslerRS2",
    "category": "method",
    "text": "Tableau for the explicit 4-stage RS2 method due to Andreas Rossler   Method cited in Table 5.3 in   Andreas Rossler, \"Second order Runge-Kutta methods for Stratonovich stochastic differential equations\",   BIT Numerical Mathematics (2007) 47   According to the paper, the method has weak order 2.0.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauSRK3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSRK3",
    "category": "method",
    "text": "Gauss-Legendre Runge-Kutta, s=3\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauSRKw1",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSRKw1",
    "category": "function",
    "text": "Tableau for the 1-stage SRKw1 method due to Wang, Hong & Xu   Method cited in   Wang, Hong, Xu, \"Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems\",   Commun. Comput. Phys. 21(1), 2017   According to the paper, the method has weak order 1.0.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauSRKw2",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSRKw2",
    "category": "function",
    "text": "Tableau for the 4-stage SRKw2 method due to Wang, Hong & Xu   Method cited in   Wang, Hong, Xu, \"Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems\",   Commun. Comput. Phys. 21(1), 2017   According to the paper, the method has weak order 2.0 when applied to systems   driven by one-dimensional noise.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauStochasticDIRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauStochasticDIRK",
    "category": "function",
    "text": "Tableau for the 2-stage stochastic symplectic DIRK method   Tableau for the stochastic symplectic DIRK method   Satisfies the conditions for Lagrange-d\'Alembert integrators.   Satisfies the conditions for strong convergence of order 1.0 for one Wiener process\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauStochasticGLRK-Tuple{Int64}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauStochasticGLRK",
    "category": "method",
    "text": "Tableau for the s-stage Gauss-Lobatto SFIRK method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauStochasticHeun-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauStochasticHeun",
    "category": "method",
    "text": "Tableau for the explicit 2-stage stochastic Heun method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauStochasticLobIIIABD2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauStochasticLobIIIABD2",
    "category": "method",
    "text": "Tableau for the 2-stage stochastic LobattoIIIA-IIIB-IIID method   Tableau for the 2-stage stochastic LobattoIIIA-IIIB-IIID method   (based on the deterministic LobattoIIIA-IIIB-IIID due to L. Jay)   It satisfies the conditions for convergence of order 1.0 for one Wiener process,   but it doesn\'t satisfy the conditions for Lagrange-d\'Alembert integrators\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauStochasticStormerVerlet-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauStochasticStormerVerlet",
    "category": "method",
    "text": "Tableau for the 2-stage stochastic LobattoIIA-IIB method (Stormer-Verlet)\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauStochasticSymplecticEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauStochasticSymplecticEuler",
    "category": "method",
    "text": "Tableau for the stochastic symplectic Euler method   Tableau for the stochastic symplectic Euler method   Satisfies the conditions for Lagrange-d\'Alembert integrators.   Satisfies the conditions for strong convergence of order 1.0 for one Wiener process   for special choices of the stochastic Hamiltonians and forces, e.g., h=h(q), f=0.\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauSymplecticEulerA-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSymplecticEulerA",
    "category": "method",
    "text": "Tableau for symplectic Euler-A method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauSymplecticEulerB-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSymplecticEulerB",
    "category": "method",
    "text": "Tableau for symplectic Euler-B method\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPGLRK-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPGLRK",
    "category": "method",
    "text": "Tableau for variational Gauss-Legendre method with s stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIA2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIA3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIA4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIA4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIA2",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIA method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIA3",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIA method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIA4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIA4",
    "category": "method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIA method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIB2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIB2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIB3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIB3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIB4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIB4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIC2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIC2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIC-III method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIC3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIC3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIC-III method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIC4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIC4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIC-III method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIID2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIID2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIID method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIID3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIID3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIID method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIID4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIID4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIID method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIE2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIE2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIE method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIE3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIE3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIE method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIE4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIE4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIE method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIF2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIF2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIF method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIF3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIF3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIF method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIF4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIF4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIF method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIG2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIG2",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIG method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIG3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIG3",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIG method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPLobIIIG4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIG4",
    "category": "method",
    "text": "Tableau for variational Gauss-Lobatto IIIG method with four stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPRadIIAIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPRadIIAIIA2",
    "category": "method",
    "text": "Tableau for Gauss-Radau IIA-IIA method with two stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPRadIIAIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPRadIIAIIA3",
    "category": "method",
    "text": "Tableau for Gauss-Radau IIA-IIA method with three stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#GeometricIntegrators.Tableaus.getTableauVPSRK3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPSRK3",
    "category": "method",
    "text": "Tableau for variational symmetric Runge-Kutta method with 3 stages\n\n\n\n\n\n"
},

{
    "location": "modules/tableaus/#Tableaus-1",
    "page": "Tableaus",
    "title": "Tableaus",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Tableaus]\nOrder   = [:constant, :type, :macro, :function]"
},

]}

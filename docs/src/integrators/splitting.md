# Splitting and Composition Methods

GeometricIntegrators supports splitting and composition methods, where the solution to an ODE of the form
```math
\dot{x} = v_1 (t,x) + ... + v_r (t,x)
```
is obtained by consecutively integrating each vector field $v_i$ independently and combining the resulting solutions in an appropriate way.
Consider a simple ODE $\dot{x} = V$ where the vector field $v$ can be written as $V = \sum_i v_i$.
The flow (exact solution) of this ODE is
```math
x(t) = \phi_{t} (x(0)) = \exp(t V) (x(0)) ,
```
and the composition method 
```math
\varphi_{\tau} = \exp(\tau v_1) \exp(\tau v_2) \dotsc \exp(\tau v_r) ,
```
where $\tau$ denotes the time step, provides a first-order accurate approximation to the exact flow as
```math
\phi_{\tau} = \exp \bigg( \tau \sum_i v_i \bigg) + \mathcal{O} (\tau^2) .
```

In the following, we use "splitting methods" to denote integrators that utilize the exact solution of each vector field $v_i$ and "composition methods" to denote integrators that utilize some consistent but possibly approximate solution for each of the vector fields $v_i$, i.e., that solution can be exact or obtained by some other integrator.
For reference see the excellent review paper by [McLachlanQuispel:2002](@cite) or the canonical book on Geometric Numerical Integration by [HairerLubichWanner:2006](@cite).


## Splitting Integrators

In GeometricIntegrators, basic splitting methods are implemented in [`IntegratorSplitting`](@ref),
which has two constructors:
```julia
IntegratorSplitting{DT,D}(solutions::Tuple, f::Vector{Int}, c::Vector, Δt)
IntegratorSplitting(equation::SODE, tableau::AbstractTableauSplitting, Δt)
```
In the first constructor, `DT` is the data type of the state vector and `D`
the dimension of the system. In the second constructor, this information
is extracted from the `equation`. 
The tuple `solutions` contains functions implementing the flow (exact solution)
of the vector fields `v_i`. The vectors `f` and `c` define the actual splitting
method: `f` is a vector of indices of the flows in the split equation to be
solved and `c` is a vector of the same size `f` that contains the coefficients
for each splitting step, i.e., the resulting integrator has the form
```math
\varphi_{\tau} = \phi_{c[s] \tau}^{v_{f[s]}} \circ \dotsc \circ \phi_{c[2] \tau}^{v_{f[2]}} \circ \phi_{c[1] \tau}^{v_{f[1]}} .
```
In the second constructor, these vectors are constructed from the tableau and
the equation.


## Composition Integrators

Fully flexible composition methods are implemented in [`IntegratorComposition`](@ref),
which can use any ODE integrator implemented in GeometricIntegrators to solve the
steps of the splitting. For each step, a different integrator can be chosen as
well as the exact solution using [`IntegratorExactODE`](@ref), which is a simple
wrapper around the exact flow of a splitting step, implementing the general
integrator interface.

`IntegratorComposition` has three constructors:
```julia
IntegratorComposition{DT,D}(integrators::Tuple, Δt)
IntegratorComposition(equation::SODE, constructors::Tuple, tableau::AbstractTableauSplitting, Δt)
IntegratorComposition(equation::SODE, tableau::AbstractTableauSplitting, Δt)
```
In the first constructor, `DT` is the data type of the state vector and `D`
the dimension of the system. In the second and third constructor, this
information is extracted from the equation. 
The tuple `integrators` contains the integrators for each substep. Each integrator
is instantiated with appropriately scaled time step size $\Delta t = c_i \tau$ to
match the corresponding splitting scheme.
In the second constructor, the tuple `constructors` contains closures around the
constructors for the integrators of each step of the composition, that is functions
taking a vector field $v_i$, the time step $\Delta t$ and optional keyword arguments,
e.g. for the exact solution or a Runge-Kutta integrator, we have
```julia
(v::Function, Δt::Number; kwargs...) -> IntegratorExactODE{DT,D}(v, Δt; kwargs...)
(v::Function, Δt::Number; kwargs...) -> Integrator{DT,D}(v, tableau, Δt; kwargs...)
```
The integrators are constructed according to the tableau and time step `\Delta t`
and passed to the first constructor.
The third constructor assumes that the exact solution is used for each splitting
step. It thus constructs a composition method that is equivalent to a plain
[`IntegratorSplitting`](@ref).


## Splitting of Hamiltonian Systems

For Hamiltonian systems, splitting is a simple and versatile technique for the construction of symplectic integrators.
Suppose that the Hamiltonian $H$ can be split into the sum of $r \geq 2$ Hamiltonians $H_{i}$ with $1 \leq i \leq r$, i.e.,
```math
H (z) = \sum \limits_{i=1}^{r} H_{i} (z) ,
```
with each Hamiltonian vector field
```math
\dot{z} = \Omega^{-T} \, \nabla H_{i} (z) 
```
explicitly solvable.
The exact solution $\phi_{t}^{H_{i}}$ of each subsystem provides a symplectic map. As the composition of symplectic maps yields a symplectic map, a symplectic integrator can be obtained by an appropriate composition of the flow maps of each subsystem.
A first-order symplectic integrator is obtained from the Lie-Trotter splitting,
```math
\varphi_{\tau} = \phi_{\tau}^{H_{1}} \circ \phi_{\tau}^{H_{2}} \circ \dotsc \circ \phi_{\tau}^{H_{k}} .
```
Second-order symplectic integrators are obtained from symmetric splittings,
```math
\varphi_{\tau} = \phi_{h/2}^{H_{1}} \circ \phi_{h/2}^{H_{2}} \circ \dotsc \circ \phi_{h/2}^{H_{k-1}} \circ \phi_{\tau}^{H_{k}} \circ \phi_{h/2}^{H_{k-1}} \circ \dotsc \circ \phi_{h/2}^{H_{2}} \circ \phi_{h/2}^{H_{1}} ,
```
or
```math
\varphi_{\tau} = \phi_{h/2}^{H_{k}} \circ \phi_{h/2}^{H_{k-1}} \circ \dotsc \circ \phi_{h/2}^{H_{2}} \circ \phi_{\tau}^{H_{1}} \circ \phi_{h/2}^{H_{2}} \circ \dotsc \circ \phi_{h/2}^{H_{k-1}} \circ \phi_{h/2}^{H_{k}} .
```
Higher order integrators can be constructed by using the Baker-Campbell-Hausdorff formula.


## Separable Hamiltonian Systems

If we have a Hamiltonian of the form $H(p,q) = T(p) + U(q)$, we can consider only the subsystem with Hamiltonian $U(q)$,
```math
\begin{aligned}
\dot{q} &= 0 , &
\dot{p} &= - \nabla U(q) .
\end{aligned}
```
This system can be solved exactly. The exact flow is
```math
\phi^{U}_{t} (q,p) = \begin{pmatrix}
q \\
p - t \nabla U(q)
\end{pmatrix} .
```
Next, consider the subsystem with Hamiltonian $T(p) = \tfrac{1}{2} p^T M^{-1} p$, 
```math
\begin{aligned}
\dot{q} &= M^{-1} p , &
\dot{p} &= 0 . 
\end{aligned}
```
This system can be solved exactly as well. The exact flow is
```math
\phi^{T}_{t} (q,p) = \begin{pmatrix}
q + t M^{-1} p \\
p
\end{pmatrix} .
```
As $\phi^{U}_{t}$ and $\phi^{T}_{t}$ are exact flows of the respective Hamiltonian, they are both symplectic.
We see that the compositions of $\phi^{U}_{t}$ and $\phi^{T}_{t}$ correspond to the symplectic Euler methods for separable Hamiltonians,
```math
\begin{aligned}
\varphi^{A}_{\tau} &= \phi^{U}_{\tau} \circ \phi^{T}_{\tau} , &
\varphi^{B}_{\tau} &= \phi^{T}_{\tau} \circ \phi^{U}_{\tau} ,
\end{aligned}
```
where $\varphi^{A}_{\tau}$ and $\varphi^{B}_{\tau}$ denote the numerical flows of symplectic Euler-A and symplectic Euler-B, respectively.
As the Störmer-Verlet methods are compositions of the symplectic Euler methods, they are also splitting methods, corresponding to the compositions
```math
\begin{aligned}
\varphi_{\tau}^{SV1} &= \varphi^{A}_{h/2} \circ \varphi^{B}_{h/2} = \phi^{U}_{h/2} \circ \phi^{T}_{\tau} \circ \phi^{U}_{h/2} , \\
\varphi_{\tau}^{SV2} &= \varphi^{B}_{h/2} \circ \varphi^{A}_{h/2} = \phi^{T}_{h/2} \circ \phi^{U}_{\tau} \circ \phi^{T}_{h/2} ,
\end{aligned}
```
respectively.
This particular splitting is often referred to as *Strang splitting* [[Strang:1968](@cite), see also [Marchuk:1968](@cite)].

Let us note that not all symplectic integrators can be obtained as splitting methods. For the symplectic Euler methods and the Störmer-Verlet methods, this is only possible for separable Hamiltonian systems. For general Hamiltonians, these methods cannot be obtained from any splitting but are nevertheless symplectic.


## Fourth Order Methods

The general form of a fourth order symplectic integrator for separable Hamiltonian systems is given by
```math
\begin{aligned}
q_{1} &= q_{0} + b_{1} \tau \, T_{p} (p_{0}) , &
p_{1} &= p_{0} - \hat{b}_{1} \tau \, U_{q} (q_{1}) , \\
q_{2} &= q_{1} + b_{2} \tau \, T_{p} (p_{1}) , &
p_{2} &= p_{1} - \hat{b}_{2} \tau \, U_{q} (q_{2}) , \\
q_{3} &= q_{2} + b_{3} \tau \, T_{p} (p_{2}) , &
p_{3} &= p_{2} - \hat{b}_{3} \tau \, U_{q} (q_{3}) , \\
q_{4} &= q_{3} + b_{4} \tau \, T_{p} (p_{3}) , &
p_{4} &= p_{3} - \hat{b}_{4} \tau \, U_{q} (q_{4}) .
\end{aligned}
```
The quantities $(q_{0}, p_{0})$ are initial values and $(q_{4}, p_{4})$ are the numerical solution after one time step $\tau$.
The whole algorithm can be written as
```math
\begin{aligned}
\varphi_{\tau} &= 
\varphi_{\hat{b}_{4} \tau}^{U}
\circ
\varphi_{b_{4} \tau}^{T}
\circ
\varphi_{\hat{b}_{3} \tau}^{U}
\circ
\varphi_{b_{3} \tau}^{T}
\circ
\varphi_{\hat{b}_{2} \tau}^{U}
\circ
\varphi_{b_{2} \tau}^{T}
\circ
\varphi_{\hat{b}_{1} \tau}^{U}
\circ
\varphi_{b_{1} \tau}^{T}
\end{aligned}
```
and is therefore immediately seen to be symplectic.
Two methods of fourth order are given by
```math
\begin{aligned}
b_{1} &= b_{4} = \dfrac{1}{2 (2 - \gamma)} , &
b_{2} &= b_{3} = \dfrac{1-\gamma}{2 (2 - \gamma)} , \\
\hat{b}_{1} &= \hat{b}_{3} = \dfrac{1}{2 - \gamma} , &
\hat{b}_{2} &= - \dfrac{\gamma}{2 - \gamma} , &
\hat{b}_{4} &= 0 ,
\end{aligned}
```
and
```math
\begin{aligned}
b_{1} &= 0 , &
b_{2} &= b_{4} = \dfrac{1}{2 - \gamma} , &
b_{3} &= \dfrac{1}{1 - \gamma^{2}} , \\
\hat{b}_{1} &= \hat{b}_{4} = \tfrac{1}{6} (2 + \gamma + \gamma^{-1}) , &
\hat{b}_{2} &= \hat{b}_{3} = \tfrac{1}{6} (2 - \gamma - \gamma^{-1}) ,
\end{aligned}
```
where $\gamma = 2^{1/3}$.
Both methods are explicit and symmetric as either $\varphi_{\hat{b}_{4} \tau}^{U}$ or $\varphi_{b_{1} \tau}^{T}$ corresponds to the identity.


## Higher Order Methods by Composition

The composition of a one-step symplectic integrator $\varphi_{\tau}$ with different step sizes provides a simple way of obtaining higher order schemes.
We assume that the initial scheme $\varphi_{\tau}$ is symmetric, that is
```math
\varphi_{\tau} \circ \varphi_{-\tau} = \mathrm{id} ,
```
as this simplifies the construction.
A symmetric method can always be built by combining a non-symmetric method with its adjoint.
If a numerical method $\varphi_{\tau}$ is symmetric, it can be used to compose higher order methods by splitting up each timestep into $s$ substeps [[HairerLubichWanner:2006](@cite), [McLachlan:1995](@cite), [MarsdenWest:2001](@cite)],
```math
\varphi_{\tau} = \varphi_{\gamma_{s} \tau} \circ ... \circ \varphi_{\gamma_{i} \tau} \circ ... \circ \varphi_{\gamma_{1} \tau} ,
```
where the careful selection of the $\gamma_{i}$ is crucial for the performance of the resulting scheme.
In the following, we present some fourth and sixth order composition methods that can be applied in most situations.


#### Fourth Order Composition Methods

If $\varphi_{\tau}$ is a method of order $r$, a method $\hat{\varphi}_{\tau}$ of order $r+2$ is obtained by the composition [[HairerLubichWanner:2006](@cite), Section V.3.2]
```math
\begin{aligned}
\hat{\varphi}_{\tau} &= \varphi_{\gamma \tau} \circ \varphi_{(1-2\gamma) \tau} \circ \varphi_{\gamma \tau} &
& \text{with} &
\gamma &= (2 - 2^{1/(r+1)})^{-1} . &
\end{aligned}
```
Hence, if $\varphi_{\tau}$ is of second order, the resulting method $\hat{\varphi}_{\tau}$ will be of fourth order.
Note that symmetric methods are always of even order.
A method of the same order but with generally smaller errors is obtained by considering five steps
```math
\begin{aligned}
\hat{\varphi}_{\tau} &= \varphi_{\gamma \tau} \circ \varphi_{\gamma \tau} \circ \varphi_{(1-4\gamma) \tau} \circ \varphi_{\gamma \tau} \circ \varphi_{\gamma \tau} &
& \text{with} &
\gamma &= (4 - 4^{1/(r+1)})^{-1} . &
\end{aligned}
```
Multiple application of these compositions yields methods of orders higher than four.


#### Sixth Order Composition Methods

Higher order compositions can also be constructed directly [[HairerLubichWanner:2006](@cite), Section V.3.2].
A sixth order method with seven substeps is given by
```math
\begin{aligned}
\gamma_{1} = \gamma_{7} &= + 0.78451361047755726381949763 , \\
\gamma_{2} = \gamma_{6} &= + 0.23557321335935813368479318 , \\
\gamma_{3} = \gamma_{5} &= - 1.17767998417887100694641568 , \\
\gamma_{4} &= + 1.31518632068391121888424973 ,
\end{aligned}
```
but again smaller errors can be achieved by using nine steps
```math
\begin{aligned}
\gamma_{1} = \gamma_{9} &= + 0.39216144400731413927925056 , \\
\gamma_{2} = \gamma_{8} &= + 0.33259913678935943859974864 , \\
\gamma_{3} = \gamma_{7} &= - 0.70624617255763935980996482 , \\
\gamma_{4} = \gamma_{6} &= + 0.08221359629355080023149045 , \\
\gamma_{5} &= + 0.79854399093482996339895035 .
\end{aligned}
```
The computational effort of these high order methods is quite large. Each step requires the solution of a nonlinear system of equations.
Given the outstanding performance already second order symplectic integrators are able to deliver, the necessity for such high order methods is rarely found.
Nevertheless, if extremely high accuracy is indispensable, think for example of long-time simulations of the solar system, these methods can be applied. Moreover, there exist special methods optimized for such problems.


## Splitting Tableaus

Actualy splitting methods are usually prescribed in one of the following forms.


### [`TableauSplitting`](@ref)

Tableau for general splitting methods for vector fields with two terms
$v = v_A + v_B$, leading to the following integrator:
```math
\varphi_{\tau} = \varphi_{b_s \tau}^{B} \circ \varphi_{a_s \tau}^{A} \circ \dotsc \circ \varphi_{b_1 \tau}^{B} \circ \varphi_{a_1 \tau}^{A} .
```


### [`TableauSplittingNS`](@ref)

Tableau for non-symmetric splitting methods [[McLachlanQuispel:2002](@cite), Equation (4.10)].
Here, two flows $\varphi_{\tau}^{A}$ and $\varphi_{\tau}^{B}$ are constructed as the Lie
composition of all vector fields in the SODE and its adjoint, respectively, i.e..
```math
\begin{aligned}
\varphi_{\tau}^{A} &= \varphi_{\tau}^{v_1} \circ \varphi_{\tau}^{v_2} \circ \dotsc \circ \varphi_{\tau}^{v_{r-1}} \circ \varphi_{\tau}^{v_r} , \\
\varphi_{\tau}^{B} &= \varphi_{\tau}^{v_r} \circ \varphi_{\tau}^{v_{r-1}} \circ \dotsc \circ \varphi_{\tau}^{v_2} \circ \varphi_{\tau}^{v_1} ,
\end{aligned}
```
and the integrator is composed as follows:
```math
\varphi_{\tau}^{NS} = \varphi_{b_s \tau}^{B} \circ \varphi_{a_s \tau}^{A} \circ \dotsc \circ \varphi_{b_1 \tau}^{B} \circ \varphi_{a_1 \tau}^{A} .
```


### [`TableauSplittingGS`](@ref)

Tableau for symmetric splitting methods with general stages [[McLachlanQuispel:2002](@cite), Equation (4.11)],
where again two flows $\varphi_{\tau}^{A}$ and $\varphi_{\tau}^{B}$ are constructed via Lie composition
```math
\begin{aligned}
\varphi_{\tau}^{A} &= \varphi_{\tau}^{v_1} \circ \varphi_{\tau}^{v_2} \circ \dotsc \circ \varphi_{\tau}^{v_{r-1}} \circ \varphi_{\tau}^{v_r} , \\
\varphi_{\tau}^{B} &= \varphi_{\tau}^{v_r} \circ \varphi_{\tau}^{v_{r-1}} \circ \dotsc \circ \varphi_{\tau}^{v_2} \circ \varphi_{\tau}^{v_1} ,
\end{aligned}
```
but with an integrator composed as
```math
\varphi_{\tau}^{GS} = \varphi_{a_1 \tau}^{A} \circ \varphi_{b_1 \tau}^{B} \circ \dotsc \circ \varphi_{b_1 \tau}^{B} \circ \varphi_{a_1 \tau}^{A} .
```


### [`TableauSplittingSS`](@ref)

Tableau for symmetric splitting methods with symmetric stages [[McLachlanQuispel:2002](@cite), Equation (4.6)].
Here, only one flow $\varphi_{\tau}^{S}$ is constructed via symmetric Strang composition,
```math
\varphi_{\tau}^{S} = \varphi_{\tau/2}^{v_1} \circ \varphi_{\tau/2}^{v_2} \circ \dotsc \circ \varphi_{\tau/2}^{v_{r-1}} \circ \varphi_{\tau/2}^{v_r} \circ \varphi_{\tau/2}^{v_r} \circ \varphi_{\tau/2}^{v_{r-1}} \circ \dotsc \circ \varphi_{\tau/2}^{v_2} \circ \varphi_{\tau/2}^{v_1} ,
```
and composed as
```math
\varphi_{\tau}^{SS} = \varphi_{a_1 \tau}^{S} \circ \varphi_{a_2 \tau}^{S} \circ \dotsc \circ \varphi_{a_s \tau}^{S} \circ \dotsc \circ \varphi_{a_2 \tau}^{S} \circ \varphi_{a_1 \tau}^{S} ,
```
to obtain an integrator.


### Implemented Splitting Methods

| Function                        | Order | Method                       |
|:--------------------------------|:------|:-----------------------------|
| [`LieA`](@ref)                  | 1     | Lie-Trotter Splitting A      |
| [`LieB`](@ref)                  | 1     | Lie-Trotter Splitting B      |
| [`Strang`](@ref)                | 2     | Strang / Marchuk Splitting   |
| [`Marchuk`](@ref)               | 2     | Strang / Marchuk Splitting   |
| [`StrangA`](@ref)               | 2     | Strang / Marchuk Splitting A |
| [`StrangB`](@ref)               | 2     | Strang / Marchuk Splitting B |
| [`McLachlan2`](@ref)            | 2     | McLachlan's 2nd order symmetric, minimum error composition method |
| [`McLachlan4`](@ref)            | 2     | McLachlan's 4th order symmetric, minimum error composition method |
| [`TripleJump`](@ref)            | 4     | 4th order "Triple Jump" composition method                        |
| [`SuzukiFractal`](@ref)         | 4     | Suzuki's 4th order "fractal" composition method                   |

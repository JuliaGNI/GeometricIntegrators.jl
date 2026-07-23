```@meta
CurrentModule = GeometricIntegrators.Integrators
```

# [Degenerate Variational Integrators](@id dvi)

Many important systems in physics are described by *degenerate* Lagrangians, in
particular Lagrangians that are linear in the velocities,
```math
\begin{equation}\label{eq:dvi-degenerate-lagrangian}
L (q, \dot{q}) = \vartheta (q) \cdot \dot{q} - H(q) ,
\end{equation}
```
where $\vartheta$ is in general a nonlinear function of $q$. The Hessian
$\partial^2 L / \partial \dot{q} \, \partial \dot{q}$ vanishes identically, so
the Lagrangian is maximally degenerate. Such Lagrangians arise, for example, for
Lotka–Volterra models, various nonlinear oscillators, and guiding-centre and
other reduced charged-particle dynamics.

Standard variational (VPRK) integrators, see the [VPRK page](@ref vprk), are
applicable when $\vartheta$ is a linear function of $q$, but are often found to be
unstable when $\vartheta$ is nonlinear [[Kraus:2017](@cite)]. Two remedies exist:
projected variational integrators (see the projection methods), and the
*degenerate variational integrators* described here, which are constructed
directly for Lagrangians of the form \eqref{eq:dvi-degenerate-lagrangian}.

GeometricIntegrators.jl provides the following degenerate variational
integrators:

| Method              | Description                                                | Order |
|:--------------------|:-----------------------------------------------------------|:------|
| [`DVIA`](@ref)      | Symplectic-Euler-A type degenerate variational integrator  | 1     |
| [`DVIB`](@ref)      | Symplectic-Euler-B type degenerate variational integrator  | 1     |
| [`CMDVI`](@ref)     | Midpoint (centred) degenerate variational integrator       | 2     |
| [`CTDVI`](@ref)     | Trapezoidal degenerate variational integrator              | 2     |
| [`DVRK`](@ref)      | Degenerate variational Runge–Kutta method                  | 2s    |


## Euler-Type and Midpoint/Trapezoidal Integrators

For definiteness we split the coordinates as $q = (x^1, x^2)$, where only the
first half of the components appear in the symplectic potential, so that the
Lagrangian takes the form
```math
L (x^1, x^2, \dot{x}^1, \dot{x}^2) = \vartheta_1 (x^1, x^2) \cdot \dot{x}^1 - H(x^1, x^2) .
```
The degenerate variational integrators are obtained from a discrete action
principle that preserves this degeneracy [[Ellison:2018](@cite)]. Applying the
discrete (Type I) phase-space action principle and introducing the velocity
$v^1_n = (x^1_{n+1} - x^1_{n}) / h$ leads to the equations of motion
```math
\begin{equation}\label{eq:dvi-eqs-of-motion}
\begin{aligned}
\dfrac{\vartheta_1 (x^1_{n}, x^2_{n}) - \vartheta_1 (x^1_{n-1}, x^2_{n-1})}{h}
  &= \nabla_{1} \vartheta_1 (x^1_{n}, x^2_{n}) \cdot v^1_{n} - \nabla_{1} H (x^1_{n}, x^2_{n}) , \\
v^1_{n} &= \big( \nabla_{2} \vartheta_1 (x^1_{n}, x^2_{n}) \big)^{-1} \nabla_{2} H (x^1_{n}, x^2_{n}) , \\
\dfrac{x^1_{n+1} - x^1_{n}}{h} &= v^1_{n} ,
\end{aligned}
\end{equation}
```
under an appropriate invertibility condition on $\nabla_2 \vartheta_1$. The
[`DVIA`](@ref) and [`DVIB`](@ref) integrators correspond to the two
symplectic-Euler-type discretisations of this principle (evaluating $\vartheta$
and $H$ at the beginning or the end of the interval) and are first order. The
[`CMDVI`](@ref) and [`CTDVI`](@ref) integrators use, respectively, a midpoint and
a trapezoidal discretisation and are symmetric and second order.

These integrators are only applicable to Lagrangians of the special form above
and are of low order; a variational construction of higher-order methods for
general degenerate Lagrangians is not currently known [[Ellison:2018](@cite)].
A higher-order, non-variational alternative is provided by the degenerate
variational Runge–Kutta methods described next.


## Degenerate Variational Runge–Kutta Methods

The [`DVRK`](@ref) methods are symplectic Runge–Kutta integrators applicable to
the degenerate Lagrangians \eqref{eq:dvi-degenerate-lagrangian}
[[Kraus:2019](@cite)]. Given a Runge–Kutta tableau with coefficients $a_{ij}$,
weights $b_{i}$ and nodes $c_{i}$, the internal stages are
```math
\begin{equation}\label{eq:dvrk}
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) = \vartheta (Q_{n,i}) , \\
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, F_{n,j} , &
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) ,
\end{aligned}
\end{equation}
```
and the update reads
```math
\begin{aligned}
q^{\mu}_{n+1} &= q^{\mu}_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V^{\mu}_{n,i} , &
p^{\mu}_{n+1} &= p^{\mu}_{n} + h \sum \limits_{i=1}^{s} b_{i} \, F^{\mu}_{n,i} ,
&& \mu = 1, \, ..., \, d/2 , \\
&&
p^{\mu}_{n+1} &= \vartheta^{\mu} (q_{n+1}) , && \mu = 1, \, ..., \, d .
\end{aligned}
```
For systems of the form \eqref{eq:dvi-degenerate-lagrangian} this method is
symplectic provided the coefficient matrix $A = (a_{ij})$ is invertible, the
coefficients satisfy the symplecticity relation
```math
\begin{equation}\label{eq:dvrk-symplecticity}
b_{i} a_{ij} + b_{j} a_{ji} = b_{i} b_{j}
\qquad \text{for all} \qquad i, j = 1, \, ..., \, s ,
\end{equation}
```
and the momentum is initialised consistently with $p_{0} = \vartheta (q_{0})$
[[Kraus:2019](@cite)]. Using Gauss–Legendre coefficients (`DVRK(Gauss(s))`)
yields a method of order $2s$.

!!! note
    The single-tableau form \eqref{eq:dvrk}–\eqref{eq:dvrk-symplecticity} follows
    [[Kraus:2019](@cite)]. The implementation accepts any
    [`Tableau`](@ref) and additionally supports a second (embedded) set of
    coefficients $\bar{a}_{ij}, \bar{b}_{i}$, for which the general symplecticity
    conditions read $b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} = b_{i} \bar{b}_{j}$
    and $\bar{b}_{i} = b_{i}$; these reduce to \eqref{eq:dvrk-symplecticity} when
    $\bar{a} = a$ and $\bar{b} = b$.

A `DVRK` integrator is constructed by passing a Runge–Kutta tableau or method:
```
DVRK(tableau::Tableau)
DVRK(method::RKMethod)
```

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
| [`DVRK`](@ref)      | Degenerate variational Runge–Kutta method                  | 2s¹   |

¹ With Gauss–Legendre coefficients, and provided the Lagrangian is given in a
gauge in which $d/2$ components of $\vartheta$ vanish identically, i.e. in the
form \eqref{eq:dvi-split-lagrangian}. See *Order of accuracy and the role of the
gauge* below.


## Euler-Type and Midpoint/Trapezoidal Integrators

For definiteness we split the coordinates as $q = (x^1, x^2)$, where only the
first half of the components appear in the symplectic potential, so that the
Lagrangian takes the form
```math
\begin{equation}\label{eq:dvi-split-lagrangian}
L (x^1, x^2, \dot{x}^1, \dot{x}^2) = \vartheta_1 (x^1, x^2) \cdot \dot{x}^1 - H(x^1, x^2) ,
\end{equation}
```
that is, $\vartheta_{\mu} = 0$ for $\mu > d/2$. Note that this is a genuine
restriction beyond \eqref{eq:dvi-degenerate-lagrangian}: $\vartheta$ is only
defined up to an exact one-form, and whether a given system is written in such a
gauge is a property of the formulation, not of the dynamics.

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
the degenerate Lagrangians \eqref{eq:dvi-split-lagrangian}
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
Note that the quadrature updates apply only to the first $d/2$ components. The
remaining components of $q_{n+1}$ receive no quadrature update and are determined
implicitly by the constraint $p_{n+1} = \vartheta (q_{n+1})$, which is also what
defines $p_{n+1}$ on output.

Symplecticity of the method rests on the coefficients satisfying the relation
```math
\begin{equation}\label{eq:dvrk-symplecticity}
b_{i} a_{ij} + b_{j} a_{ji} = b_{i} b_{j}
\qquad \text{for all} \qquad i, j = 1, \, ..., \, s ,
\end{equation}
```
together with three hypotheses on the problem and one on the initial data. In full,
the method preserves the noncanonical symplectic two-form
$\omega = \mathrm{d}\vartheta$ provided that [[Kraus:2019](@cite)]

1. the Lagrangian is given in the form \eqref{eq:dvi-split-lagrangian}, i.e. in a
   gauge with $\vartheta_{\mu} = 0$ for $\mu > d/2$;
2. the $d/2 \times d/2$ matrix $\partial \vartheta_{\mu} / \partial q^{\nu}$, with
   $\mu \le d/2 < \nu$, is invertible;
3. the coefficient matrix $A = (a_{ij})$ is invertible;
4. the coefficients satisfy \eqref{eq:dvrk-symplecticity};
5. the momentum is initialised consistently with $p_{0} = \vartheta (q_{0})$.

Condition 2 is what makes $\omega$ nondegenerate — under condition 1 one has
$\det \omega = (\det \partial \vartheta_{\mu} / \partial q^{\nu})^{2}$ — and what
makes the scheme well posed, since it is exactly what allows
$p_{n+1} = \vartheta (q_{n+1})$ to be solved for the components of $q_{n+1}$ that
carry no quadrature update.

!!! warning
    Conditions 1 and 2 are independent, and condition 2 is the weaker of the two.
    For `LotkaVolterra2d`, whose gauge $\vartheta = (q_{2} + \log q_{2} / q_{1}, \,
    q_{1})$ violates condition 1, the matrix of condition 2 is
    $\partial \vartheta_{1} / \partial q^{2} = 1 + 1 / q_{1} q_{2} \ne 0$ and hence
    perfectly invertible. Invertibility therefore does *not* certify that the
    problem is in class; see *Order of accuracy and the role of the gauge* below.

The `DVRK` constructor checks conditions 3 and 4 and warns if either is violated
(pass `check_conditions = false` to suppress the warnings); condition 5 is checked
when the integrator cache is built. Conditions 1 and 2 are properties of the
problem and are not checked.

!!! note "Admissible tableaus"
    Conditions \eqref{eq:dvrk-symplecticity} and the invertibility of $A$ are
    jointly restrictive. They exclude Lobatto IIIA (singular $A$) and Radau IIA
    (invertible $A$, but \eqref{eq:dvrk-symplecticity} is violated), and they
    cannot be met by the Lobatto IIIA–IIIB pair, which is symplectic only as a
    *partitioned* method. The Gauss–Legendre tableaus satisfy both for any number
    of stages, and are the only family used with `DVRK` in practice.

### Order of accuracy and the role of the gauge

`DVRK(Gauss(s))` attains the full order $2s$ — but only for Lagrangians that
genuinely satisfy \eqref{eq:dvi-split-lagrangian}, i.e. for which $d/2$
components of $\vartheta$ vanish identically. The symplectic potential $\vartheta$
of a given system is only defined up to an exact one-form, and the *gauge* matters:

| Problem | Gauge | $s = 1$ | $s = 2$ | $s = 3$ |
|:--------|:------|:-------:|:-------:|:-------:|
| `LotkaVolterra2dSingular`         | $\vartheta_{2} = 0$ (in class)   | 2.00 | 4.00 | 5.98 |
| `LotkaVolterra2d`                 | $\vartheta_{2} = q_{1} \ne 0$    | 1.01 | 2.00 | 3.01 |
| `LotkaVolterra2dSymmetric`        | both components $\ne 0$          | 0.99 | 2.00 | 3.01 |
| `MasslessChargedParticleSingular` | $A_{2} = 0$ (in class)           | 2.00 | 4.00 | 5.98 |
| `MasslessChargedParticle`         | both components $\ne 0$          | 1.00 | 2.00 | 3.00 |

The rows within each block describe the *same dynamics*, differing only by a gauge
transformation of $\vartheta$. Applied outside the class, `DVRK` remains
convergent, but its order is halved from $2s$ to $s$. When using `DVRK`, choose the
formulation of the problem in which the appropriate components of $\vartheta$
vanish identically — the `…Singular` variants in `GeometricProblems` exist for
exactly this purpose.

!!! warning
    Deleting a component of $\vartheta$ is not a gauge transformation. For the
    Lotka–Volterra model, $\vartheta = (q_{2} + \log q_{2} / q_{1}, \, q_{1})$ and
    $\vartheta = (\log q_{2} / q_{1}, \, 0)$ differ by the exact one-form
    $\mathrm{d}(q_{1} q_{2})$ and describe the same system, whereas
    $(q_{2} + \log q_{2} / q_{1}, \, 0)$ describes a *different* system — one whose
    orbits coincide but whose time parametrisation does not.

A `DVRK` integrator is constructed by passing a Runge–Kutta tableau or method:
```
DVRK(tableau::Tableau; check_conditions = true)
DVRK(method::RKMethod; check_conditions = true)
```

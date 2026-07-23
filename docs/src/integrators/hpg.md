```@meta
CurrentModule = GeometricIntegrators.Integrators
```

# [Hamilton-Pontryagin-Galerkin Integrators](@id hpg)

!!! note
    The Hamilton–Pontryagin–Galerkin (HPG) framework described on this page is
    not yet implemented in GeometricIntegrators.jl. This page documents the
    theory that underlies the closely related, finite-difference
    [Hamilton–Pontryagin integrators](@ref hpi) (`HPImidpoint`, `HPItrapezoidal`),
    which arise as particular low-order instances of the framework.

Hamilton–Pontryagin–Galerkin integrators are obtained from a discontinuous
Galerkin discretisation of the Hamilton–Pontryagin principle. Whereas the
construction of a variational integrator usually involves two ingredients — a
finite-dimensional function space that approximates the phase-space trajectories,
and a quadrature rule that approximates the action integral — the HPG framework
adds a third ingredient, a *continuity constraint* or *jump condition* that
connects the (possibly broken) function spaces across time nodes. Many known
variational integrators, obtained from distinct variational principles, emerge as
special cases of this unifying framework, including the discrete Euler–Lagrange
equations, the position-momentum form, the symplectic Euler and Störmer–Verlet
methods, and symplectic and variational Runge–Kutta methods.


## Continuous Hamilton–Pontryagin Principle

The Hamilton–Pontryagin principle [[Yoshimura:2006b](@cite)] is an action
principle on the Pontryagin bundle $TQ \oplus T^*Q$. In contrast to Hamilton's
principle of stationary action, where the variations $\delta v$ are induced by
variations $\delta q$, here the variations of the velocity $v$ are left free.
Instead, a kinematic constraint enforcing the second-order condition $v = \dot{q}$
is added, with the momentum $p$ acting as a Lagrange multiplier,
```math
\begin{equation}\label{eq:hpg-continuous-action}
\mathcal{A} [q,v,p] = \int \limits_{0}^{T} \Big[ L(q,v) + \left< p , \dot{q} - v \right> \Big] \, dt .
\end{equation}
```
Computing variations,
```math
\delta \mathcal{A} [q,v,p] = \int \limits_{0}^{T} \bigg[ \frac{\partial L}{\partial q} \, \delta q + \frac{\partial L}{\partial v} \, \delta v + \left< \delta p , \dot{q} - v \right> + \left< p , \delta \dot{q} - \delta v \right> \bigg] \, dt = 0 ,
```
yields the second-order condition, the Legendre transform, and the
Euler–Lagrange equations,
```math
\begin{aligned}
\dot{q} &= v , &
p &= \frac{\partial L}{\partial v} (q,v) , &
\dot{p} &= \frac{\partial L}{\partial q} (q,v) .
\end{aligned}
```
Equivalently, introducing the generalised energy $E(q,v,p) = \left< p, v \right> - L(q,v)$, the principle can be written as
```math
\mathcal{A} [q,v,p] = \int \limits_{0}^{T} \Big[ \left< p, \dot{q} \right> - E(q,v,p) \Big] \, dt ,
```
whose variations lead to a set of generalised Hamilton equations
```math
\begin{aligned}
\dot{q} &= \frac{\partial E}{\partial p} (q,v,p) , &
\frac{\partial E}{\partial v} &= 0 , &
\dot{p} &= - \frac{\partial E}{\partial q} (q,v,p) .
\end{aligned}
```
For constrained systems, the Hamilton–Pontryagin principle generalises to the
Lagrange–d'Alembert–Pontryagin principle [[Yoshimura:2006a](@cite), [Yoshimura:2006b](@cite)].


## Galerkin Discretisation

The trajectories $\big( q(t), v(t), p(t) \big)$ on the Pontryagin bundle are
approximated by piecewise polynomials on a fixed time interval $\mathcal{I} = [0, T]$,
partitioned into $N$ subintervals $\mathcal{I}_{n} = (t_{n}, t_{n+1})$ with
$t_{n} = nh$ and fixed time step $h$. On each subinterval, $q$, $v$ and $p$ are
represented by polynomials of degree $r_q$, $r_v$ and $r_p$,
```math
Q_{n} (t) = \sum \limits_{i=1}^{r_{q}} Q_{n,i} \, \varphi_{n,i}^{q} (t) ,
\qquad
V_{n} (t) = \sum \limits_{i=1}^{r_{v}} V_{n,i} \, \varphi_{n,i}^{v} (t) ,
\qquad
P_{n} (t) = \sum \limits_{i=1}^{r_{p}} P_{n,i} \, \varphi_{n,i}^{p} (t) ,
```
with basis functions $\varphi_{n,i}$, e.g. Lagrange polynomials, and quadrature
nodes $c_{i}$ with weights $b_{i}$ (satisfying $\sum_{i=1}^{s} b_{i} = 1$),
located at $t_{n,i} = t_{n} + h c_{i}$. With this data the discrete
Hamilton–Pontryagin action \eqref{eq:hpg-continuous-action} becomes
```math
\begin{multline}\label{eq:hpg-discrete-action}
\mathcal{A}_{d} [Q,V,P] = h \sum \limits_{n=0}^{N-1} \Bigg[
\sum \limits_{i=1}^{s} b_{i} \Big( L \big( Q_{n} (t_{n,i}) , V_{n} (t_{n,i}) \big)
+ \big< P_{n} (t_{n,i}) , \, \dot{Q}_{n} (t_{n,i}) - V_{n} (t_{n,i}) \big> \Big) \\
+ \text{(continuity constraints or numerical flux)} \Bigg] .
\end{multline}
```
Since the function spaces may be broken, the values of the polynomials $Q_n$ and
$Q_{n+1}$ at a shared node $t_{n+1}$ need not coincide; they are denoted by the
one-sided limits $q_{n+1}^{-}$ and $q_{n+1}^{+}$, and analogously for $v$ and $p$.
For continuous discretisations $q_{n+1} = q_{n+1}^{-} = q_{n+1}^{+}$, whereas for
discontinuous discretisations one typically chooses the average
$q_{n+1} = \tfrac{1}{2} (q_{n+1}^{-} + q_{n+1}^{+})$.


## Continuity Constraints and Numerical Fluxes

The last term in \eqref{eq:hpg-discrete-action} encodes how neighbouring
subintervals are connected. Different choices reproduce the different types of
generating functions of classical mechanics:

| Constraint | Continuity of $p$    | Continuity of $q$    |
|:-----------|:---------------------|:---------------------|
| $(q,Q)$    | double discontinuous | left-right-continuous|
| $(q,P)$    | right-continuous     | left-continuous      |
| $(p,Q)$    | left-continuous      | right-continuous     |
| $(p,P)$    | left-right-continuous| double discontinuous |

Equivalently, the constraints can be written in terms of a numerical flux
$\hat{p}_{n+1} \, \hat{q}_{n+1}$ with parameters $\alpha, \beta \in [0,1]$,
```math
\begin{aligned}
\hat{p}_{n+1} &= (1 - \alpha) \, P_{n} (t_{n+1}) + \alpha \, P_{n+1} (t_{n+1}) , &
\hat{q}_{n+1} &= (1 - \beta) \, Q_{n} (t_{n+1}) + \beta \, Q_{n+1} (t_{n+1}) ,
\end{aligned}
```
which interpolate between the $(q,P)$ constraint ($\alpha = 0$, $\beta = 1$) and
the $(p,Q)$ constraint ($\alpha = 1$, $\beta = 0$). By varying the polynomial
degrees, quadrature rule and continuity constraint, this framework recovers many
established variational integrators and also yields new classes of high-order,
discontinuous-Galerkin integrators. The finite-difference midpoint and
trapezoidal Hamilton–Pontryagin integrators implemented in this package (see the
[Hamilton–Pontryagin page](@ref hpi)) correspond to the lowest-order members of
this family; Hamilton–Pontryagin integrators on Lie groups are discussed in
[[BouRabee:2009](@cite)], and the connection to discrete Dirac mechanics in
[[Leok:2011](@cite)].

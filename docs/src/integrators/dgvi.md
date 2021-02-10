# [Discontinuous Galerkin Variational Integrators](@id dgvi)

Discontinuous Galerkin Variational Integrators (DGVIs) are a family of integrators
for degenerate Lagrangian systems and for Hamiltonian systems subject to Dirac constraints.
For integrators for non-degenerate (regular) Lagrangian and unconstrained
Hamiltonian systems see [Hamilton-Pontryagin-Galerkin (HPG) Integrators](@ref hpg).


## Degenerate Lagrangian Systems

Consider a fully degenerate Lagrangian system of the form
```math
L(q, \dot{q}) = \vartheta (q) \cdot \dot{q} - H(q) ,
```
where $\vartheta (q)$ denotes the Cartan one-form and $H(q)$ the Hamiltonian,
that is usually given by the total energy of the system.


## Discrete Trajectories and Numerical Quadrature

The first step in the derivation of variational integrators is the discretization
of the action integral
```math
\mathcal{A} [q] = \int \limits_{0}^{T} L(q(t), \dot{q}(t)) \, dt .
```
To this end, the interval $[0,T]$ is split into $N$ sub-intervals
$[t_{n}, t_{n+1}]$ with $t_n = nh$ and $h$ the time step size, so that $t_N = T$
and the action can be written as
```math
\mathcal{A} [q] = \sum \limits_{n=0}^{N-1} \int \limits_{[t_{n}, t_{n+1}]} L(q(t), \dot{q}(t)) \, dt .
```
Within each interval $(t_{n}, t_{n+1})$ a piecewise-polynomial approximation $q_h$
of the trajectory $q$ is constructed using $S$ basis functions $\varphi_{i}$,
```math
q_h(t) \vert_{(t_{n}, t_{n+1})} = \sum \limits_{i=1}^{S} x_{n,i} \, \bar{\varphi}_{n,i} (t) ,
```
where $\bar{\varphi}_{n,i} (t)$ is a rescaled basis function, defined by
```math
\bar{\varphi}_{n,i} (t) = \varphi_{i} \bigg( \frac{t - t_{n}}{t_{n+1} - t_{n}} \bigg) ,
```
and it is assumed that $\varphi_{i}$ is compactly supported on $[0,1]$.
These approximations $q_h(t)$ do not need to be continuous across interval boundaries
but are indeed allowed to have jumps.
Replacing the continuous trajectory $q$ in the action with $q_h$, we obtain
```math
\mathcal{A} [q_h] = \sum \limits_{n=0}^{N-1} \int \limits_{(t_{n}, t_{n+1})} \big[ \vartheta (q_h (t)) \cdot \dot{q}_h (t) - H(q_h (t)) \big] \, dt
+ \sum \limits_{n=1}^{N-1} [[ \vartheta (q_h (t)) \cdot \dot{q}_h (t) ]]_{t=t_n}  .
```
The integral of the Hamiltonian $H(q_h)$ over the interval boundaries does not
contribute to the integral, differently from the term $\vartheta (q_h) \cdot \dot{q_h}$,
which will determine the numerical flux $[[ \cdot ]]_n$ at $t_n$ of the Discontinuous Galerkin method.
The approximation of this term will be discussed below.
In order to obtain a fully discrete action, a numerical quadrature rule with $R$
nodes $c_i$ and weights $b_i$ is introduced for the approximation of the integral,
```math
\mathcal{A}_d [x_d] = h \sum \limits_{n=0}^{N-1} \sum \limits_{i=1}^{R} b_i \big[ \vartheta (q_h(t_n + c_i h)) \cdot \dot{q}_h (t_n + c_i h) - H(q_h(t_n + c_i h)) \big]
+ \sum \limits_{n=1}^{N-1} [[ \vartheta (q_h (t)) \cdot \dot{q}_h (t) ]]_{t=t_n} .
```
Here, $x_d$ denotes the vector of all the degrees of freedom, i.e.,
```math
x_d = ( x_{0,1}, ..., x_{0,S}, \, x_{1,1}, ..., x_{N-2,S}, \, x_{N-1,1}, ..., x_{N-1,S} )^T .
```
In order to write the discrete action in a more explicit form, mass and derivative
matrices $m$ and $a$ are introduced, whose elements are given by
```math
m_{ij} = \varphi_j (c_i) ,
\qquad
a_{ij} = \varphi_j' (c_i) ,
\qquad
i = 1, ..., R ,
\;
j = 1, ..., S ,
```
so that the solution and its time derivative at the quadrature points can be written as
```math
Q_{n,i} \equiv q_h(t_n + c_i h) = m_{ij} x_{n,j} ,
\qquad
V_{n,i} \equiv \dot{q}_h (t_n + c_i h) = a_{ij} x_{n,j} ,
```
where
```math
x_{n} = ( x_{n,1}, ..., x_{n,S} )^T
```
is the vector containing the degrees of freedom of $q_h \vert_{[t_{n}, t_{n+1}]}$.
Using these definitions, the discrete action can be written as
```math
\mathcal{A}_d [x_d] = h \sum \limits_{n=0}^{N-1} \sum \limits_{i=1}^{R} b_i \big[ \vartheta (Q_{n,i}) \cdot V_{n,i} - H(Q_{n,i}) \big]
+ \sum \limits_{n=1}^{N-1} [[ \vartheta (q_h (t)) \cdot \dot{q}_h (t) ]]_{t=t_n} .
```


## Numerical Fluxes

In the following, the solution values "left" and "right" of the jump will be needed.
This will be denoted by $q_n^-$ and $q_n^+$, respectively.
Usually, these just correspond to the polynomials on the left and right, evaluated at $t_n$, i.e.,
```math
q_n^- = \lim_{t \uparrow t_n} q_h (t) = q_h \vert_{[t_{n-1}, t_{n}]} (t_n) ,
\qquad
q_n^+ = \lim_{t \downarrow t_n} q_h (t) = q_h \vert_{[t_{n}, t_{n+1}]} (t_n) .
```
In principle, however, more general reconstructions of the solution could be used.
In the following, it will be assumed, that $q_n^-$ is given by some linear
combinations of the degrees of freedom of the polynomial on the left interval
and correspondingly that $q_n^+$ is given by some linear
combinations of the degrees of freedom of the polynomial on the right interval,
specifically
```math
q_n^- = r^{-} \cdot x_{n} ,
\qquad
q_n^+ = r^{+} \cdot x_{n+1} ,
```
where $r^{\pm}$ are appropriate coefficient vectors.


#### Gauge Terms

The Lagrangian $L$ can be augmented with any total time derivative without
changing the (continuous) Euler-Lagrange equations. In particular, one can
consider the modified Lagrangian
```math
L'(q, \dot{q}) = \vartheta (q) \cdot \dot{q} - H(q) - \nu \dfrac{d}{dt} \bigg( \vartheta (q) \cdot q \bigg) .
```
While this gauge term vanishes in the continuous case, it takes a finite
value across jumps of the discontinuous discrete solution, so that the modified
discrete action reads
```math
\mathcal{A}_d' [x_d] = h \sum \limits_{n=0}^{N-1} \sum \limits_{i=1}^{R} b_i \big[ \vartheta (Q_{n,i}) \cdot V_{n,i} - H(Q_{n,i}) \big]
+ \sum \limits_{n=1}^{N-1} \bigg[\bigg[ \vartheta (q_h (t)) \cdot \dot{q}_h (t) - \nu \dfrac{d}{dt} \bigg( \vartheta (q) \cdot q \bigg) \bigg]\bigg]_{t=t_n} .
```
In the following, only the modified Lagrangian and action will be considered,
in order to obtain a sufficiently general framework for constructing numerical
fluxes.
For brevity of notation, the prime will be dropped.


#### Total Time Derivatives Across Jumps

The computation of the total time derivative in the gauge term is simple, at
least in the distributional sense.
Even though both, $\vartheta (q_h)$ and $q_h$ have a jump, the jump occurs at the
same position in time, so that the derivative can be computed as
```math
\dfrac{d}{dt} \bigg( \vartheta (q_h) \cdot q_h \bigg) \bigg\vert_{t=t_n}
= \dfrac{d}{dt} \bigg( \vartheta (q_n^-) \cdot q_n^- \, \Theta (t_n - t) + \vartheta (q_n^+) \cdot q_n^+ \, \Theta (t - t_n) \bigg) ,
```
where $\Theta$ denotes the Heaviside function.
This can be explicitly computed as
```math
\dfrac{d}{dt} \bigg( \vartheta (q_h) \cdot q_h \bigg) \bigg\vert_{t=t_n}
= - \vartheta (q_n^-) \cdot q_n^- \, \delta (t_n) + \vartheta (q_n^+) \cdot q_n^+ \, \delta (t_n)
= [[ \vartheta (q_h) \cdot q_h ]]_{t=t_n} \, \delta (t_n) ,
```
with $\delta (t_n)$ the Dirac delta-function at $t_n$.


#### Non-conservative Products

Simple means for integrating the Lagrangian across jumps are provided by
discretisations of the integral
```math
\int \limits_0^1 \vartheta (\Phi(\tau; q^-, q^+)) \cdot \dfrac{d \Phi(\tau; q^-, q^+)}{d\tau} \, d\tau ,
```
where $\Phi$ is a path connecting the solution values $q^-$ and $q^+$ on the
left and the right of the jump.
Upon picking a quadrature rule with $\sigma$ nodes $\gamma_i$ and corresponding
weights $\beta_i$, the discrete product takes the form
```math
\sum \limits_{i=1}^{\sigma} \beta_i \, \vartheta \big( \Phi (\gamma_i; q_{n}^- , q_{n}^+) \big) \cdot \dfrac{d\Phi}{d\tau} (\gamma_i; q_{n}^- , q_{n}^+) .
```
For a compact notation, "mass" and "derivative" vectors $\mu^{\pm}$ and $\alpha^{\pm}$ are introduced, so that
```math
\Phi (\gamma_i; q_{n}^- , q_{n}^+) = \mu^-_i q_{n}^- + \mu^+_i q_{n}^+
\qquad
\Phi' (\gamma_i; q_{n}^- , q_{n}^+) = \alpha^-_i q_{n}^- + \alpha^+_i q_{n}^+ ,
```
and the discrete product can be written as
```math
\sum \limits_{i=1}^{\sigma} \beta_i \, \vartheta ( \mu^-_i q_{n}^- + \mu^+_i q_{n}^+ ) \cdot ( \alpha^-_i q_{n}^- + \alpha^+_i q_{n}^+ ) .
```
Providing the path $\Phi$ by two functions $\phi^{\pm}(\tau)$, so that
```math
\phi(\tau; q^-, q^+) = q^- \phi^{-}(\tau) + q^+ \phi^{+}(\tau) ,
```
the components of the "mass" and "derivative" vectors are given by
```math
\mu^-_i = \phi^{-} (\gamma_i) ,
\qquad
\mu^+_i = \phi^{+} (\gamma_i) ,
```
and
```math
\alpha^-_i = \frac{d\phi^{-}}{d\tau} (\gamma_i) ,
\qquad
\alpha^+_i = \frac{d\phi^{+}}{d\tau} (\gamma_i) ,
```
respectively.


## Discrete Variational Principle

Using the construction of the previous sections, the discrete action reads
```math
\mathcal{A}_d [x_d] = h \sum \limits_{n=0}^{N-1} \sum \limits_{i=1}^{R} b_i \big[ \vartheta (Q_{n,i}) \cdot V_{n,i} - H(Q_{n,i}) \big] \\
+ \sum \limits_{n=1}^{N-1} \bigg[ \sum \limits_{i=1}^{\sigma} \beta_i \, \vartheta ( \mu^-_i q_{n}^- + \mu^+_i q_{n}^+ ) \cdot ( \alpha^-_i q_{n}^- + \alpha^+_i q_{n}^+ )
- \nu \big[ \vartheta (q_n^+) \cdot q_n^+ - \vartheta (q_n^-) \cdot q_n^- \big] \bigg] .
```
The discrete Euler-Lagrange equations are obtained by applying Hamilton's
principle of stationary action to $\mathcal{A}_d [x_d]$, that is requiring that
$\delta \mathcal{A}_d [x_d] = 0$.
The variations of the discrete action are computed as follows,
```math
\delta \mathcal{A}_d [x_d]
= h \sum \limits_{n=0}^{N-1} \sum \limits_{i=1}^{R} b_i \big[ \delta Q_{n,i} \cdot \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} + \vartheta (Q_{n,i}) \cdot \delta V_{n,i} - \delta Q_{n,i} \cdot \nabla H(Q_{n,i}) \big] \\
+ \sum \limits_{n=1}^{N-1} \bigg[ \sum \limits_{i=1}^{\sigma} \beta_i \, \big[ ( \mu^-_i \delta q_{n}^- + \mu^+_i \delta q_{n}^+ ) \cdot \nabla \vartheta ( \mu^-_i q_{n}^- + \mu^+_i q_{n}^+ ) \cdot ( \alpha^-_i q_{n}^- + \alpha^+_i q_{n}^+ ) + \vartheta ( \mu^-_i q_{n}^- + \mu^+_i q_{n}^+ ) \cdot ( \alpha^-_i \delta q_{n}^- + \alpha^+_i \delta q_{n}^+ ) \big] \\
- \nu \big[ \delta q_n^+ \cdot \nabla \vartheta (q_n^+) \cdot q_n^+  + \vartheta (q_n^+) \cdot \delta q_n^+ - \delta q_n^- \cdot \nabla \vartheta (q_n^-) \cdot q_n^- - \vartheta (q_n^-) \cdot \delta q_n^- \big] \bigg] .
```
Using the relations
```math
\delta Q_{n,i} = m_{ij} \delta x_{n,j} ,
\qquad
\delta V_{n,i} = \frac{a_{ij}}{h} \delta x_{n,j} ,
\qquad
\delta q_{n}^- = r^{-}_{j} \delta x_{n-1,j} ,
\qquad
\delta q_{n}^+ = r^{+}_{j} \delta x_{n,j} ,
```
the variations of the discrete action become
```math
\delta \mathcal{A}_d [x_d]
= \sum \limits_{n=0}^{N-1} \sum \limits_{i=1}^{R} \sum \limits_{j=1}^{S} b_i \big[ h m_{ij} \delta x_{n,j} \cdot \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} + \vartheta (Q_{n,i}) \cdot a_{ij} \delta x_{n,j} - h m_{ij} \delta x_{n,j} \cdot \nabla H(Q_{n,i}) \big] \\
+ \sum \limits_{n=1}^{N-1} \bigg[ \sum \limits_{i=1}^{\sigma} \sum \limits_{j=1}^{S} \beta_i \, \big[ ( \mu^-_i r^{-}_{j} \delta x_{n-1,j} + \mu^+_i r^{+}_{j} \delta x_{n,j} ) \cdot \nabla \vartheta ( \mu^-_i q_{n}^- + \mu^+_i q_{n}^+ ) \cdot ( \alpha^-_i q_{n}^- + \alpha^+_i q_{n}^+ ) + \vartheta ( \mu^-_i q_{n}^- + \mu^+_i q_{n}^+ ) \cdot ( \alpha^-_i r^{-}_{j} \delta x_{n-1,j} + \alpha^+_i r^{+}_{j} \delta x_{n,j} ) \big] \\
- \nu \big[ r^{+}_{j} \delta x_{n,j} \cdot \nabla \vartheta (q_n^+) \cdot q_n^+  + \vartheta (q_n^+) \cdot r^{+}_{j} \delta x_{n,j} - r^{-}_{j} \delta x_{n-1,j} \cdot \nabla \vartheta (q_n^-) \cdot q_n^- - \vartheta (q_n^-) \cdot r^{-}_{j} \delta x_{n-1,j} \big] \bigg] .
```
Requiring the variation of the discrete action to vanish yields the discrete equations of motion,
```math
0 = \sum \limits_{i=1}^{R} b_i \big[ h m_{ij} \nabla \vartheta (Q_{n,i}) \cdot V_{n,i} + a_{ij} \vartheta (Q_{n,i}) - h m_{ij} \nabla H(Q_{n,i}) \big] \\
+ \bigg[ \sum \limits_{i=1}^{\sigma} \beta_i \, \big[
   \mu^-_i r^{-}_{j} \nabla \vartheta ( \mu^-_i q_{n+1}^- + \mu^+_i q_{n+1}^+ ) \cdot ( \alpha^-_i q_{n+1}^- + \alpha^+_i q_{n+1}^+ )
 + \mu^+_i r^{+}_{j} \nabla \vartheta ( \mu^-_i q_{n  }^- + \mu^+_i q_{n  }^+ ) \cdot ( \alpha^-_i q_{n  }^- + \alpha^+_i q_{n  }^+ ) \\
 + \vartheta ( \mu^-_i q_{n+1}^- + \mu^+_i q_{n+1}^+ ) \alpha^-_i r^{-}_{j}
 + \vartheta ( \mu^-_i q_{n  }^- + \mu^+_i q_{n  }^+ ) \alpha^+_i r^{+}_{j}
   \big] \\
 - \nu \big[
     r^{+}_{j} \nabla \vartheta (q_{n  }^+) \cdot q_{n  }^+ + r^{+}_{j} \vartheta (q_{n  }^+)
   - r^{-}_{j} \nabla \vartheta (q_{n+1}^-) \cdot q_{n+1}^- - r^{-}_{j} \vartheta (q_{n+1}^-)
   \big]
 \bigg] ,
```
for all $n$ and all $j$.

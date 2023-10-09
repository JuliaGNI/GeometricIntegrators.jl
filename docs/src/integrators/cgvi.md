# [Continuous Galerkin Variational Integrators](@id cgvi)

In the following, we first present the Galerkin framework for variational integrators [[Leok:2004](@cite), [Leok:2011](@cite), [Campos:2014](@cite), [Campos:2015](@cite), [OberBloebaum:2014](@cite), [OberBloebaum:2016](@cite)] where the space of curves $\mf{Q}$ that connect two points in $\mf{M}$ is approximated by a finite-dimensional subspace (Galerkin integrators of 0th kind).
We use Lagrange polynomials to approximate the trajectories albeit other choices are possible.
Alternatively, we can approximate the generalised velocities 
(Galerkin integrators of 1st kind), which leads us to variational-partitioned Runge-Kutta methods.
Even though, in some cases correspondences between integrators of 0th kind and integrators of 1st kind might be found, in general they are different.

In our treatment we distinguish between integer timesteps, which are the coordinates $q_{n}$ at time $t_{n}$, and internal stages (nodes), which are the coordinates $Q_{n,i}$ located between two integer timesteps $q_{n}$ and $q_{n+1}$ at consecutive points in time $t_{n}$ and $t_{n+1}$.
In all of the following we assume that the timestep $h$ is constant.

In the discrete variational principle, we have to consider variations at both, the integer timesteps and the internal stages.
For Galerkin integrators of 0th kind, the coordinates $Q_{n,i}$ are varied, whereas for the Galerkin integrators of 1st kind, the velocities $\dot{Q}_{n,i}$ are varied. So the independent variables are $(q_{n}, Q_{n,i})$ and $(q_{n}, \dot{Q}_{n,i})$, respectively.


## Space of Discrete Trajectories

In order to construct the discrete space of curves from $q_{0}$ to $q_{N}$,
```math
\mf{Q}_{d} ( q_{0}, q_{N}, \{ t_{n} \}_{n=0}^{N} ) = \big\{ q_{d} : \{ t_{n} \}_{n=0}^{N} \rightarrow \mf{M} \; \big\vert \; q_{d} (t_{0}) = q_{0}, \, q_{d} (t_{N}) = q_{N} \big\} \subset \mf{Q}_{d} ( \mf{M} ) ,
```
we will construct a finite-dimensional subspace of
```math
\mf{Q} ( q_{0}, q_{N}, [0, T] ) = \big\{ q : [0, T] \rightarrow \mf{M} \; \big\vert \; q(0) = q_{0}, \, q(T) = q_{N} \big\} \subset \mf{Q} ( \mf{M} ) .
```
The subspace $\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] )$ is defined by
```math
\begin{equation}\label{eq:galerkin-vi-space-of-piecewise-polynomials}
\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] ) = \big\{q_{h} : [0, T] \rightarrow \mf{M} \; \big\vert \; q_{h} \vert_{ [ t_{n} , t_{n+1} ] } \in \mathbb{P}_{s} ([t_{n}, t_{n+1}]) , \;  q_{h} \in C^{0} ([0, T]) \big\} ,
\end{equation}
```
where $\mathbb{P}_{s} ([t_{n}, t_{n+1}])$ is the space of polynomials of degree $s$ in the interval $[t_{n}, t_{n+1}] \subset [0, T]$.
We see that $\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] ) \subset \mf{Q} ( q_{0}, q_{N}, [0, T] )$. In order for $\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] )$ to be an instance of $\mf{Q}_{d} ( q_{0}, q_{N}, \{ t_{n} \}_{n=0}^{N} )$ we have to require in addition to the definition \eqref{eq:galerkin-vi-space-of-piecewise-polynomials} that  on the sequence $\{ t_{n} \}_{n=0}^{N}$, the curves $q_{h}$ in $\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] )$ satisfy
```math
q_{h} \vert_{ [ t_{n} , t_{n+1} ] } (t_{n})   = q_{n}
\qquad \text{and} \qquad
q_{h} \vert_{ [ t_{n} , t_{n+1} ] } (t_{n+1}) = q_{n+1} ,
```
where $q_{n}$ are the points of the discrete trajectories $q_{d} = \{ q_{n} \}_{n=0}^{N}$.
However, it is often more convenient to enforce these continuity constraints weakly in the action rather than building them into the spaces, which implies dropping the condition $q_{h} \in C^{0} ([0, T])$ in \eqref{eq:galerkin-vi-space-of-piecewise-polynomials}. This in turn means that the whole of $\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] )$ is not a subspace of $\mf{Q} ( q_{0}, q_{N}, [0, T] )$ anymore.

**TODO**:
*The role of continuity is not clearly explained and probably not correctly worked out.
In practice, we are using broken spaces $\mf{Q}_h$ and enforce continuity only weakly in the
action principle. 
It should be made clear, that $\mf{Q}_d$ is an approximation of $\mf{Q}$, but not a subspace of $\mf{Q}$.
$\mf{Q}_h$, however, is a subspace of $\mf{Q}$ (although it can be broken).
The connection between $\mf{Q}_h$ and $\mf{Q}_d$ is made by the continuity constraints (which for
Lagrange polynomials and sequences of nodes which include the boundaries, i.e., $c_1=0$
and $c_s=1$, is automatically satisfied -> this is not true).
In this context we also need to discuss that curves in $\mf{Q}(\mf{M})$ are assumed to be $C^2$,
which is not required by the curves in $\mf{Q}_h$.*


## Piecewise Polynomials

A basis of $\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] )$ can easily be constructed by combining bases of $\mathbb{P}_{s} ([t_{n}, t_{n+1}])$, e.g., using Lagrange polynomials.
We start by specifying the collocation times of the internal stages. Select a set of $s$ points $c_{i}$ with $0 \leq c_{i} \leq 1$, which are the nodes of the basis functions. Therefore the basis is built by $s$ functions.
The internal stages are then located at $t_{n,i} = t_{n} + h c_{i}$, such that
```math
t_{n} \leq t_{n,1} < ... < t_{n,s} \leq t_{n+1} ,
```
and we have
```math
Q_{n,i} \approx q ( t_{n,i} ) .
```
We will only be concerned with Lagrange polynomials.
The $j$-th Lagrange polynomial of order $s$ is defined by
```math
l^{s,i} (\tau) = \prod \limits_{\substack{1 \leq j \leq s,\\ j \neq i}} \dfrac{\tau - c_{j}}{c_{i} - c_{j}} .
```
The $c_{i}$ are often chosen to be the collocation points of some quadrature rule (e.g., Gauß-Legendre or Gauß-Chebyshew points).
Within each subinterval between two consecutive timesteps, $t_{n}$ and $t_{n+1}$, the same Lagrange basis is used, namely
```math
\spn \big\{ \varphi_{n}^{s,m} (t) \; \big\vert \; 1 \leq m \leq s \big\} ,
```
with
```math
\varphi_{n}^{s,m} (t) = \begin{cases}
l^{s,m} \big( (t-t_{n}) / (t_{n+1} - t_{n}) \big) , \hspace{5em} & t_{n} \leq t \leq t_{n+1} , \\
\hspace{.8em} 0 , & \text{else} .
\end{cases}
```
It suffices to specify the basis for one subinterval $[t_{n}, t_{n+1}]$ and then replicate this basis for all subintervals, so that the finite-dimensional subspace of $\mf{Q} ( q_{0}, q_{N}, [0, T] )$ is given by
```math
\begin{multline}
\mf{Q}_{h} ( q_{0}, q_{N}, [0, T] ) 
= \Big\{ q_{h} : [0, T] \rightarrow \mf{M} \; \Big\vert \; q_{h} \big\vert_{[t_{n}, t_{n+1}]} \in \spn \big\{ \varphi_{n}^{s,m} \; \big\vert \; 1 \leq m \leq s \big\} , \\ %\, \text{for} \, 0 \leq n \leq N-1 , \\
  q_{h} (0) = q_{0} , \, q_{h} (T) = q_{N} \Big\} .
\end{multline}
```
In order to obtain the discrete space of curves $\mf{Q}_{d} ( q_{0}, q_{N}, \{ t_{n} \}_{n=0}^{N} )$ we have to add continuity constraints, connecting the polynomials in each interval with the nodal values $q_{n}$, that is
```math
\begin{multline}
\mf{Q}_{d} ( q_{0}, q_{N}, \{ t_{n} \}_{n=0}^{N} ) 
= \Big\{ q_{h} \in \mf{Q}_{h} ( q_{0}, q_{N}, [0, T] ) \; \Big\vert \; q_{h} \big\vert_{[t_{n}, t_{n+1}]} (t_{n}) = q_{n} , \\ q_{h} \big\vert_{[t_{n}, t_{n+1}]} (t_{n+1}) = q_{n+1} \Big\} .
\end{multline}
```
We project the trajectories of the particles onto a Lagrange basis in order to obtain the polynomial approximation of the trajectory in the interval $[t_{n}, t_{n+1}]$, i.e.,
```math
\begin{equation}\label{eq:particle-trajectory-position}
q_{h} (t) \big\vert_{[t_{n}, t_{n+1}]} = \sum \limits_{m=1}^{s} Q_{n,m} \, \varphi_{n}^{s,m} (t) .
\end{equation}
```
The particle velocities are then obtained by differentiating with respect to time,
```math
\begin{equation}\label{eq:particle_trajectory_velocity}
\dot{q}_{h} (t) \big\vert_{[t_{n}, t_{n+1}]} = \sum \limits_{m=1}^{s} Q_{n,m} \, \dot{\varphi}_{n}^{s,m} (t) .
\end{equation}
```
For $s=2$, we obtain
```math
\begin{aligned}
q_{h} (t) \big\vert_{ [ t_{n} , t_{n+1} ] } &= Q_{n,1} \, \dfrac{t - t_{n+1}}{t_{n} - t_{n+1}} + Q_{n,2} \, \dfrac{t - t_{n}}{t_{n+1} - t_{n}} , &
\dot{q}_{h} (t) \big\vert_{ [ t_{n} , t_{n+1} ] } &= \dfrac{Q_{n,2} - Q_{n,1}}{t_{n+1} - t_{n}} ,
\end{aligned}
```
which is just linear interpolation for $q$ and piecewise constant for $\dot{q}$.

**TODO**:
*Discuss other basis functions (e.g. Chebychev polynomials, B-splines) and other quadrature rules (e.g., Chebyshev points, optimised IGA points) and visualise basis functions for different quadrature points.*


## Numerical Quadrature

In order to numerically compute the definite integral
```math
\begin{equation}\label{eq:galerkin-vi-quadrature-integral}
F [q] = \int \limits_{t_{n}}^{t_{n+1}} f \big( t, q(t) \big) \, dt ,
\end{equation}
```
we apply two levels of approximation. As $q(t)$ is unknown, we replace it with the piecewise polynomial approximation $q_{h} (t)$.
Further, we introduce a quadrature formula in which $f$ itself is approximated by a Lagrange polynomial, i.e.,
```math
f_{h} (t, q (t)) = \sum \limits_{i=1}^{s} \varphi_{n}^{s,i} (t) \, f \big( t_{n} + h c_{i}, \, q (t_{n} + h c_{i}) \big) .
```
Together, this gives an approximation of the integral \eqref{eq:galerkin-vi-quadrature-integral} as follows
```math
\begin{equation}\label{eq:galerkin_vi_quadrature_rule}
F_h [q]
= \int \limits_{t_{n}}^{t_{n+1}} f_{h} (t, q_{h} (t)) \, dt
= h \sum \limits_{i=1}^{s} b_{i} \, f \big( t_{n} + h c_{i}, \, q_{h} (t_{n} + h c_{i}) \big) ,
\end{equation}
```
where $b_{i}$ are the weights or coefficients of the quadrature formula, given by
```math
\begin{equation}\label{eq:galerkin_vi_quadrature_weights}
b_{i}
= \dfrac{1}{h} \int \limits_{t_{n}}^{t_{n+1}} \varphi_{n}^{s,i} (t) \, dt
= \int \limits_{0}^{1} l^{s,i} (\tau) \, d\tau ,
\end{equation}
```
$h = t_{n+1} - t_{n}$ is the time step, and $q_{h} (t) \big\vert_{ [ t_{n} , t_{n+1} ] }$ is some polynomial approximation to $q(t)$ in the interval $[t_{n}, t_{n+1}]$.
%The $c_{i}$ will also be the collocation points of the quadrature rule that is used to approximate the action integral.
We will focus on collocation methods where the nodes $c_{i}$ of the quadrature rule are also the nodes of the basis functions, so that for \eqref{eq:particle-trajectory-position}, we have $q_{h} (t_{n} + h c_{i}) = Q_{n,i}$.
It follows that the discrete Lagrangian can be written as
```math
L_{d} (q_{n}, q_{n+1}) = h \sum \limits_{i=1}^{s} b_{i} \, L \big( Q_{n,i} , \dot{Q}_{n,i} \big) .
```
Most often, we use Gauss-Legendre quadrature, where the nodes $c_{i}$ are given by the roots of the Legendre polynomials. The Gauss quadrature rules with $s$ nodes yield exact results when applied to polynomials of order up to $2s$.


## Galerkin Variational Integrators

In order to write the discrete Lagrangian in the discrete action in a compact form, we define the nodal coefficients $a_{ij}$,
```math
a_{ij}
= h \dfrac{d \varphi_{n}^{s,j}}{dt} \bigg\vert_{t=t_{n} + h c_{i}}
= \dfrac{d l^{s,j}}{d\tau} \bigg\vert_{\tau=c_{i}} ,
```
so that the velocities can be written as
```math
\dot{Q}_{n,i}
= \dot{q}_{h} (t_{n,i})
= \sum \limits_{j=1}^{s} Q_{n,j} \, \dot{\varphi}_{n}^{s,j} (t_{n} + h c_{i})
= \dfrac{1}{h} \sum \limits_{j=1}^{s} a_{ij} \, Q_{n,j} .
```
In order to complete the discrete action, we explicitly add the continuity constraint,
```math
\begin{multline}
\mathcal{A}_{d} [q_{d}]
= \sum \limits_{n=0}^{N-1} \bigg[
	\sum \limits_{i=1}^{s} b_{i} \, L \big( Q_{n,i} , \dot{Q}_{n,i} \big)
  + \lambda_{n} \cdot \big( q_{h}\vert_{ [ t_{n} , t_{n+1} ] } (t_{n  }) - q_{n  } \big) \\
  + \mu_{n+1}   \cdot \big( q_{h}\vert_{ [ t_{n} , t_{n+1} ] } (t_{n+1}) - q_{n+1} \big)
\bigg] ,
\end{multline}
```
which ensures that the polynomials in neighbouring intervals, e.g., $[t_{n}, t_{n+1}]$ and $[t_{n+1}, t_{n+2}]$, have the same value at integer timesteps, e.g., $t_{n+1}$.


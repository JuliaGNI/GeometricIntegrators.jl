```@meta
CurrentModule = GeometricIntegrators
```

# Runge-Kutta Methods

Runge-Kutta methods exploit the **Fundamental Theorem of Calculus**, which states that the solution of an initial-value problem
```math
\begin{aligned}
\dot{x} (t) &= f(t, x(t)) , &
x(t_{n}) &= x_{n} ,
\end{aligned}
```
at time $t_{n+1}$ is given by
```math
x (t_{n+1}) = x (t_{n}) + \int \limits_{t_{n}}^{t_{n+1}} \dot{x} (t) \, dt .
```

Runge-Kutta methods are constructed by approximating the integral by some quadrature formula with $s$ nodes $c_{i}$ and corresponding weights $b_{i}$ to obtain $x_{n+1} \approx x (t_{n+1})$ by
```math
\begin{aligned}
x_{n+1} &= x_{n} + h \sum \limits_{i=1}^{s} b_{i} \dot{X}_{n,i} , &
\dot{X}_{n,i} &= f(t_{n} + c_{i} h, X_{n,i}) ,
\end{aligned}
```
where the internal stage values $X_{n,i} \approx x(t_{n} + c_{i} h)$ for $i = 1, ..., s$ are determined by another quadrature formula, approximating the integral
```math
x(t_{n} + c_{i} h) = x (t_{n}) + \int \limits_{t_{n}}^{t_{n} + c_i h} \dot{x} (t) \, dt ,
```
namely
```math
X_{n,i} = x_{n} + h \sum \limits_{j=1}^{s} a_{ij} \dot{X}_{n,j} ,
```
with the same vector field values $\dot{X}_{n,j}$ used for the computation of $x_{n+1}$.

**Definition:** Runge-Kutta methods are numerical one-step methods
```math
\begin{aligned}
X_{n,i} &= x_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, f(t_{n} + c_{j} h, X_{n,j}) , \\
x_{n+1} &= x_{n} + h \sum \limits_{j=1}^{s} b_{j}  \, f(t_{n} + c_{j} h, X_{n,j}) ,
\end{aligned}
```
defined by a set of nodes $c_i$, weights $b_i$ and coefficients $a_{ij}$ with $i,j = 1, ..., s$, summarized in the Butcher tableau
```math
\begin{array}{c|c}
c & a     \\
\hline
  & b^{T} \\
\end{array}
=
\begin{array}{c|cccc}
c_{1}  & a_{11} & a_{12} & \dots & a_{1s} \\
c_{2}  & a_{21} & a_{22} & \dots & a_{2s} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
c_{s}  & a_{s1} & a_{s2} & \dots & a_{ss} \\
\hline
       & b_{1}  & b_{2}  & \dots & b_{s}  \\
\end{array}
```

Most properties of the methods, such as order or stability, can be analysed just by posing conditions on the Butcher tableau.


## Common Runge-Kutta Methods

GeometricIntegrators.jl provides various explicit and implicit (both diagonally and fully implicit) Runge-Kutta methods.
For many methods, tabulated coefficients are included, namely


| Function and Aliases                                           | Stages | Order |
|:---------------------------------------------------------------|:-------|:------|
| **Explicit Methods**                                           |        |       |
| [`TableauExplicitEuler`](@ref), [`TableauForwardEuler`](@ref)  | 1      | 1     |
| [`TableauExplicitMidpoint`](@ref)                              | 2      | 2     |
| [`TableauHeun2`](@ref)                                         | 2      | 2     |
| [`TableauHeun3`](@ref)                                         | 3      | 3     |
| [`TableauKutta`](@ref), [`TableauKutta3`](@ref)                | 3      | 3     |
| [`TableauRalston2`](@ref)                                      | 2      | 2     |
| [`TableauRalston3`](@ref)                                      | 3      | 3     |
| [`TableauRunge`](@ref), [`TableauRunge2`](@ref)                | 2      | 2     |
| [`TableauRK416`](@ref), [`TableauRK4`](@ref)                   | 4      | 4     |
| [`TableauRK438`](@ref)                                         | 4      | 4     |
| [`TableauSSPRK3`](@ref)                                        | 3      | 3     |
| **Diagonally Implicit Methods**                                |        |       |
| [`TableauCrankNicolson`](@ref)                                 | 2      | 2     |
| [`TableauCrouzeix`](@ref)                                      | 2      | 3     |
| [`TableauKraaijevangerSpijker`](@ref)                          | 2      | 2     |
| [`TableauQinZhang`](@ref)                                      | 2      | 2     |
| **Fully Implicit Methods**                                     |        |       |
| [`TableauImplicitEuler`](@ref), [`TableauBackwardEuler`](@ref) | 1      | 1     |
| [`TableauImplicitMidpoint`](@ref)                              | 2      | 2     |
| [`TableauSRK3`](@ref)                                          | 3      | 4     |

The coefficients of other methods are computed on-the-fly as described in the following.


## Simplifying Assumptions

The construction of many Runge-Kutte methods, in particular the Gauß, Lobatto and Radau methods, relies on the so-called simplifying assumptions:
```math
\begin{aligned}
B(\sigma): & \sum \limits_{i=1}^{s} b_{i} c_{i}^{k-1} = \frac{1}{k} , &
k = 1 , \, ... , \, \sigma , \\
%
C(\eta): & \sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_{i}^{k}}{k} , &
i = 1, \, ... , \, s , \; k = 1 , \, ... , \, \eta , \\
%
D(\zeta): & \sum \limits_{i=1}^{s} b_{i} c_{i}^{k-1} a_{ij} = \frac{b_{j}}{k} ( 1 - c_{j}^{k} ) , &
j = 1, \, ... , \, s , \; k = 1 , \, ... , \, \zeta .
\end{aligned}
```

These assumptions provide order conditions for a Runge-Kutta method given by $(a_{ij}, b_{i}, c_{i})$.
The condition $B(p)$ implies that the quadrature rule $(b_{i}, c_{i})$ is of order $p$.
Furthermore, the following theorem holds:

**Theorem (Butcher 1964)**: If the coefficients $(a_{ij}, b_{i}, c_{i})$ of a Runge-Kutta method satisfy $B(\sigma)$, $C(\eta)$, $D(\zeta)$ with $\sigma \le \eta + \zeta + 1$ and $\sigma \le 2 \eta + 2$, then the method is of order $\sigma$.


## Gauß, Lobatto and Radau Methods

Gauß methods are collocation methods using the nodes and weights of Gaußian quadrature formulas.
The nodes are the zeros of the shifted Legendre polynomials of degree $s$,
```math
\frac{d^{s}}{dx^{s}} \big( x^{s} (x-1)^{s} \big) .
```
In a similar fashion, the nodes of the Radau I and II and the Lobatto III methods are defined as the roots of the polynomials
```math
\begin{aligned}
& \frac{d^{s-1}}{dx^{s-1}} \big( x^s (x-1)^{s-1} \big)     && \text{(Radau I)} , \\
& \frac{d^{s-1}}{dx^{s-1}} \big( x^{s-1} (x-1)^s \big)     && \text{(Radau II)} , \\
& \frac{d^{s-2}}{dx^{s-2}} \big( x^{s-1} (x-1)^{s-1} \big) && \text{(Lobatto III)} .
\end{aligned}
```

The weights $b_{1}, ..., b_{s}$ are chosen such that the methods satisfy $B(\sigma)$, that is $B(s)$, for the Gauß methods, $B(s-1)$ for the Radau methods, and $B(s-2)$ for the Lobatto methods.

The coefficients $a_{ij}$ for $i,j = 1, ..., s$ are obtained by the simplifying assumption $C(s)$ for the Gauß, Radau IIA and Lobatto IIIA methods, and by the simplifying assumption $D(s)$ for the Radau IA and Lobatto IIIB methods.
The coefficients of the Lobatto IIIC methods are determined by setting $a_{i,1} = b_1$ for $i = 1, ..., s$ and solving the simplifying assumption $C(s-1)$, while the coefficients of the Lobatto IIIC̄ methods are determined by setting $a_{i,s} = 0$ and solving $C(s-1)$. Note that the Lobatto IIIC̄ methods are sometimes also called Lobatto III or Lobatto III*. For reasons of code symmetry we chose to stick with the less common name Lobatto IIIC̄.
The Lobatto IIID and IIIE methods are obtained by combining the tableaus of the Lobatto IIIC and IIIC̄ and the Lobatto IIIA and IIIB methods, respectively, i.e., 
```math
\begin{aligned}
a_{ij}^{D} &= \tfrac{1}{2} ( a_{ij}^{C} + a_{ij}^{C̄} ) &
& \text{and} &
a_{ij}^{E} &= \tfrac{1}{2} ( a_{ij}^{A} + a_{ij}^{B} ) .
\end{aligned}
```
While the Lobatto IIIA, IIIB, IIIC and IIIC̄ methods are not symplectic on their own (although the Lobatto IIIA-IIIB and IIIC-IIIC̄ pairs constitute symplectic partitioned Runge-Kutta methods), the Lobatto IIID and IIIE methods are each symplectic by themselves.

The Gauß methods are of order $2s$, the Radau methods or order $2s-1$ and the Lobatto methods are of order $2s-2$, with the exception of the Lobatto IIIF method. This method has been specifically constructed to be of order $2s$ as described in [[Fangzong:2016](@cite)].
The Lobatto IIIG method is constructed in a similar fashion as the Lobatto IIID and IIIE methods by averaging the coefficients of the Lobatto IIIF method with its symplectic complement, i.e.,
```math
a_{ij}^{G} = \tfrac{1}{2} ( a_{ij}^{F} + \bar{a}_{ij}^{F} ) ,
```
where the coffiecients $\bar{a}_{ij}^{F}$ are determined by
```math
\begin{aligned}
b_{i} \bar{a}^{F}_{ij} + \bar{b}_{j} a^{F}_{ji} &= b_{i} \bar{b}_{j}  &
& \text{and} &
\bar{b}_{i} &= b_{i} .
\end{aligned}
```
The tableaus of all of the above methods can be computed for an arbitrary number of stages $s$ and thus to arbitrary order.


### Constructors

The following methods are provided for selecting the previously described Runge-Kutta schemes:

| Function                         | Method                      | Order |
|:---------------------------------|:----------------------------|:------|
| [`TableauGauss(s)`](@ref)        | Gauß-Legendre with s stages | 2s    |
| [`TableauLobattoIIIA(s)`](@ref)  | Lobatto IIIA with s stages  | 2s-2  |
| [`TableauLobattoIIIB(s)`](@ref)  | Lobatto IIIB with s stages  | 2s-2  |
| [`TableauLobattoIIIC(s)`](@ref)  | Lobatto IIIC with s stages  | 2s-2  |
| [`TableauLobattoIIIC̄(s)`](@ref)  | Lobatto IIIC̄ with s stages  | 2s-2  |
| [`TableauLobattoIIID(s)`](@ref)  | Lobatto IIID with s stages  | 2s-2  |
| [`TableauLobattoIIIE(s)`](@ref)  | Lobatto IIIE with s stages  | 2s-2  |
| [`TableauLobattoIIIF(s)`](@ref)  | Lobatto IIIF with s stages  | 2s-2  |
| [`TableauLobattoIIIG(s)`](@ref)  | Lobatto IIIG with s stages  | 2s-2  |
| [`TableauRadauIA(s)`](@ref)      | Radau IA with s stages      | 2s-1  |
| [`TableauRadauIB(s)`](@ref)      | Radau IB with s stages      | 2s-1  |
| [`TableauRadauIIA(s)`](@ref)     | Radau IIA with s stages     | 2s-1  |
| [`TableauRadauIIB(s)`](@ref)     | Radau IIB with s stages     | 2s-1  |

The first argument `s` refers to the number of stages ($s \ge 1$ for Gauß and $s \ge 2$ for all other methods).
The second argument specifies the number type of the coefficients. Internally, all coefficients are computed using `BigFloat` and then converted to the requested number type, defaulting to `Float64`.


## Partitioned Equations

Partitioned Runge-Kutta methods consist of two tableaus that solve a partitioned ordinary differential equation,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{aligned}
```
in the following way:
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, v(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i}  \, v(t_{n} + c_j \Delta t, Q_{n,i}, P_{n,i}) , \\
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, f(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i}   \, f(t_{n} + c_j \Delta t, Q_{n,i}, P_{n,i}) .
\end{aligned}
```

The [`PartitionedTableau`](@ref) data structure can be used to compose any two Runge-Kutta tableaus
into a partitioned Runge-Kutta tableau.
A particular interesting family of partitioned Runge-Kutta methods are symplectic Lobatto methods,
specifically

| Function                               | Method                      | Order |
|:---------------------------------------|:----------------------------|:------|
| [`TableauLobattoIIIAIIIB(s)`](@ref)    | Lobatto-IIIA-IIIB           | 2s-2  |
| [`TableauLobattoIIIBIIIA(s)`](@ref)    | Lobatto-IIIB-IIIA           | 2s-2  |
| [`TableauLobattoIIIAIIIĀ(s)`](@ref)    | Lobatto-IIIA-IIIĀ           | 2s-2  |
| [`TableauLobattoIIIBIIIB̄(s)`](@ref)    | Lobatto-IIIB-IIIB̄           | 2s-2  |
| [`TableauLobattoIIICIIIC̄(s)`](@ref)    | Lobatto-IIIC-IIIC̄           | 2s-2  |
| [`TableauLobattoIIIC̄IIIC(s)`](@ref)    | Lobatto-IIIC̄-IIIC           | 2s-2  |
| [`TableauLobattoIIIDIIID̄(s)`](@ref)    | Lobatto-IIID-IIID̄           | 2s-2  |
| [`TableauLobattoIIIEIIIĒ(s)`](@ref)    | Lobatto-IIIE-IIIĒ           | 2s-2  |
| [`TableauLobattoIIIFIIIF̄(s)`](@ref)    | Lobatto-IIIF-IIIF̄           | 2s    |
| [`TableauLobattoIIIF̄IIIF(s)`](@ref)    | Lobatto-IIIF̄-IIIF           | 2s    |
| [`TableauLobattoIIIGIIIḠ(s)`](@ref)    | Lobatto-IIIG-IIIḠ           | 2s    |


## Implicit Equations

An implicit ordinary differential equations is an initial value problem of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) .
\end{aligned}
```
Such problems can be integrated with adapted Runge-Kutta methods, namely
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, v(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i}  \, v(t_{n} + c_j \Delta t, Q_{n,i}, P_{n,i}) , \\
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, f(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i}   \, f(t_{n} + c_j \Delta t, Q_{n,i}, P_{n,i}) , \\
P_{n,i} &= ϑ(t_{n} + c_j \Delta t, Q_{n,j}, P_{n,j}) .
\end{aligned}
```

Implicit ODEs can be integrated with any implicit Runge-Kutta or partitioned Runge-Kutta method.


## Custom Tableaus

If required, it is straight-forward to create a custom tableau.
The tableau of Heun's method, for example, is defined as follows:
```@example 1
using GeometricIntegrators # hide
a = [[0.0 0.0]
     [1.0 0.0]]
b = [0.5, 0.5]
c = [0.0, 1.0]
o = 2

tab = Tableau(:heun, o, a, b, c)
```
Here, `o` is the order of the method, `a` are the coefficients, `b` the weights
and `c` the nodes. For partitioned Runge-Kutta tableaus, `PartitionedTableau` can
be used. The first parameter of the constructor of each tableau assigns a name to
the tableau.
Such custom tableaus can be used in exactly the same as standard tableaus, making
it very easy to implement and test new Runge-Kutta methods.

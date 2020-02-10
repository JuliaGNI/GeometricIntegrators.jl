@doc raw"""
`HODE`: Hamiltonian Ordinary Differential Equation *EXPERIMENTAL*

Defines a Hamiltonian ordinary differential initial value problem, that is
a canonical Hamiltonian system of equations,
```math
\begin{align*}
\dot{q} (t) &= v(t, q(t), p(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , & p(t_{0}) &= p_{0} ,
\end{align*}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{align*}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} ,
\end{align*}
```
initial conditions ``(q_{0}, p_{0})`` and the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `n`: number of initial conditions
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `h`: function computing the Hamiltonian ``H``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``

"""

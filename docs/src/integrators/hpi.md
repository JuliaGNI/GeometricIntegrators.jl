# [Hamilton-Pontryagin Integrators](@id hpi)

Hamilton-Pontryagin integrators are obtained from discrete versions of the Hamilton-Pontryagin action principle,
```math
\delta \mathcal{A} [q,v,p] = \delta \int_{０}^{T} \left[ L (q, v) - \left< p , \dot{q} - v \right> \right] \, dt = 0 .
```

## Trapezoidal Discrete Lagrangian

We consider a discrete Lagrangian of the form
```math
L_d (q_{n}, q_{n+1}) = \frac{h}{2} \big[ L (q_{n}, v_{n+1/2}) + L (q_{n+1}, v_{n+1/2}) \big] ,
```
where $q_{n}$ approximates the solution $q(t_n)$ and $v_{n+1/2}$ is the velocity, which is assumed to be constant in the interval $[t_{n}, t_{n+1}]$.
The discrete Hamilton-Pontryagin action reads
```math
A_d [q_d] = \sum \limits_{n=0}^{N-1} \bigg[ L_d (q_{n}, q_{n+1}) + h \left< p_{n+1/2} , \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) - v_{n+1/2} \right> \bigg] ,
```
where $\phi_h$ is a map that computes the velocity $v_{n+1/2}$ as a function of $q_{n}$, $q_{n+1}$ and a set of parameters $a_{n,n+1}$.
A trivial example of such a map that does not depend on any parameters $a_{n,n+1}$ is
```math
\phi_h (q_{n}, q_{n+1}; a_{n,n+1}) = \frac{q_{n+1} - q_{n}}{h} .
```
In order to solve the discrete Euler-Lagrange equations, the user needs to specify the map $\phi_h$, its gradients with respect to $q_{n}$ and $q_{n+1}$, denoted by $D_1 \phi_h$ and $D_2 \phi_h$, respectively, the gradient with respect to the parameters, denoted by $D_a \phi_h$, and an initial set of parameters $a_{0}$.

The equations of motion, that are solved by this integrator, are computed as:
```math
\begin{aligned}
0 &= \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n}, v_{n+1/2})
   + \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n}, v_{n-1/2}) \\
  &+ h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2}
   + h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n-1/2} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
v_{n+1/2}
&= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) , \\
p_{n+1/2}
&= \frac{1}{2} \, \frac{\partial L}{\partial v} (q_{n}, v_{n+1/2})
 + \frac{1}{2} \, \frac{\partial L}{\partial v} (q_{n+1}, v_{n+1/2}) .
\end{aligned}
```
Upon defining the momentum
```math
p_{n}
= h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n-1/2}
+ \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n}, v_{n-1/2}) ,
```
we can rewrite the equations of motion as
```math
\begin{aligned}
0 &= p_{n}
   + \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n}, v_{n+1/2})
   + h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
v_{n+1/2}
&= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) , \\
p_{n+1/2}
&= \frac{1}{2} \, \frac{\partial L}{\partial v} (q_{n}, v_{n+1/2})
 + \frac{1}{2} \, \frac{\partial L}{\partial v} (q_{n+1}, v_{n+1/2}) , \\
p_{n+1}
&= h \, D_2 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2}
 + \frac{h}{2} \, \frac{\partial L}{\partial q} (q_{n+1}, v_{n+1/2}) .
\end{aligned}
```
Given $(q_{n}, p_{n})$, the first four equations can be solved for $q_{n+1}$, where $v_{n+1/2}$, $p_{n+1/2}$, and $a_{n,n+1}$ are treated as internal variables similar to the internal stages of a Runge-Kutta method, and the last equation provides an explicit update for $p_{n+1}$.


## Midpoint Discrete Lagrangian

Similarly to the previous section, we can consider a discrete Lagrangian of the form
```math
L_d (q_{n}, q_{n+1}) = h \, L \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) ,
```
where $q_{n}$ approximates the solution $q(t_n)$ and $v_{n+1/2}$ is the velocity, which is assumed to be constant in the interval $[t_{n}, t_{n+1}]$.
The discrete Hamilton-Pontryagin action reads
```math
A_d [q_d] = \sum \limits_{n=0}^{N-1} \bigg[ L_{d}^{\alpha,\beta} (q_{n}, q_{n+1}, v_{n+1/2}) + \left< p_{n+1/2} , \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) - v_{n+1/2} \right> \bigg] ,
```
where $\phi_h$ is a map that computes the velocity $v_{n+1/2}$ as a function of $q_{n}$, $q_{n+1}$ and a set of parameters $a_{n,n+1}$.

The equations of motion, that are solved by this integrator, is computed as:
```math
\begin{aligned}
0 &= \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n-1} + q_{n}}{2}, v_{n-1/2} \bigg)
   + \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) \\
  &+ h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2}
   + h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n-1/2} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
v_{n+1/2} &= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) , \\
p_{n+1/2} &= \frac{\partial L}{\partial v} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) .
\end{aligned}
```
Upon defining the momentum
```math
p_{n}
= h \, D_2 \phi_h (q_{n-1}, q_{n}; a_{n-1,n}) \cdot p_{n-1/2}
+ \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n-1} + q_{n}}{2}, v_{n-1/2} \bigg) ,
```
we can rewrite the equations of motion as
```math
\begin{aligned}
0 &= p_{n}
   + \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg)
   + h \, D_1 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
0 &= D_a \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2} , \\
v_{n+1/2} &= \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) , \\
p_{n+1/2} &= \frac{\partial L}{\partial v} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) , \\
p_{n+1}
&= h \, D_2 \phi_h (q_{n}, q_{n+1}; a_{n,n+1}) \cdot p_{n+1/2}
 + \frac{h}{2} \, \frac{\partial L}{\partial q} \bigg( \frac{q_{n} + q_{n+1}}{2}, v_{n+1/2} \bigg) .
\end{aligned}
```
Given $(q_{n}, p_{n})$, the first four equations can be solved for $q_{n+1}$, where $v_{n+1/2}$, $p_{n+1/2}$, and $a_{n,n+1}$ are treated as internal variables similar to the internal stages of a Runge-Kutta method, and the last equation provides an explicit update for $p_{n+1}$.

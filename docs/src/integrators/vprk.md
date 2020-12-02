```@meta
CurrentModule = GeometricIntegrators.Integrators.VPRK
```

# Variational Partitioned Runge-Kutta Integrators

Variational partitioned Runge-Kutta methods solve Lagranian systems in implicit form, i.e.,
```math
\begin{aligned}
p       &= \dfrac{\partial L}{\partial \dot{q}} (q, \dot{q}) , &
\dot{p} &= \dfrac{\partial L}{\partial q}       (q, \dot{q}) , 
\end{aligned}
```
by the following scheme,
```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial \dot{q}} (Q_{n,i}, V_{n,i}) , &
F_{n,i} &= \dfrac{\partial L}{\partial q}       (Q_{n,i}, V_{n,i}) , \\
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij}       \, V_{n,j} , &
P_{n,i} &= p_{n} + h \sum \limits_{j=1}^{s} \bar{a}_{ij} \, F_{n,j} , \\
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i}        \, V_{n,i} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i}  \, F_{n,i} .
\end{aligned}
```
Here, $s$ denotes the number of internal stages, $a_{ij}$ and $\bar{a}_{ij}$ are the coefficients of the Runge-Kutta method and $b_{i}$ and $\bar{b}_{i}$ the corresponding weights.
If the coefficients satisfy the symplecticity conditions,
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} &
& \text{and} &
\bar{b}_{i} &= b_{i} ,
\end{aligned}
```
these methods correspond to the position-momentum form of the discrete Lagrangian [[MarsdenWest:2001](@cite)]
```math
L_{d} (q_{n}, q_{n+1}) = h \sum \limits_{i=1}^{s} b_{i} \, L \big( Q_{n,i}, V_{n,i} \big) .
```

While these integrators show favourable properties for systems with regular Lagrangian, they are usually not applicable for degenerate Lagrangian systems, in particular those with Lagrangians of the form $L (q, \dot{q}) = \vartheta(q) \cdot \dot{q} - H(q)$.
While variational integrators are still applicable in the case of $\vartheta$ being a linear function of $q$, they are often found to be unstable when $\vartheta$ is a nonlinear function of $q$ as is the case with Lotka-Volterra systems, guiding centre dynamics and various nonlinear oscillators.
To mitigate this problem, projection methods have been developed, which when applied to variational integrators provide long-time stable integrators for general degenerate Lagrangian systems that maintain conservation of energy and momenta [[Kraus:2017](@cite)].

GeometricIntegrators.jl provides the following VPRK methods (*some are still experimental*):

| Integrator                            | Description                                                                                          |
|:--------------------------------------|:-----------------------------------------------------------------------------------------------------|
| [`IntegratorVPRK`](@ref)              | Variational Partitioned Runge-Kutta (VPRK) integrator without projection                             |
| [`IntegratorVPRKpStandard`](@ref)     | VPRK integrator with standard projection                                                             |
| [`IntegratorVPRKpSymmetric`](@ref)    | VPRK integrator with symmetric projection                                                            |
| [`IntegratorVPRKpMidpoint`](@ref)     | VPRK integrator with midpoint projection                                                             |
| [`IntegratorVPRKpVariational`](@ref)  | VPRK integrator with variational projection                                                          |
| [`IntegratorVPRKpSecondary`](@ref)    | VPRK integrator with projection on secondary constraint                                              |
| [`IntegratorVPRKpInternal`](@ref)     | Gauss-Legendre VPRK integrator with projection on internal stages of Runge-Kutta method              |
| [`IntegratorVPRKpTableau`](@ref)      | Gauss-Legendre VPRK integrator with projection in tableau of Runge-Kutta method                      |

For testing purposes [`IntegratorVPRKpStandard`](@ref) provides some additional constructors (*these methods are generally unstable*):

| Integrator                            | Description                                                                                          |
|:--------------------------------------|:-----------------------------------------------------------------------------------------------------|
| [`IntegratorVPRKpVariationalQ`](@ref) | VPRK integrator with variational projection on $(q_{n}, p_{n+1})$                                    |
| [`IntegratorVPRKpVariationalP`](@ref) | VPRK integrator with variational projection on $(p_{n}, q_{n+1})$                                    |
| [`IntegratorVPRKpSymplectic`](@ref)   | VPRK integrator with symplectic projection                                                           |

All of the above integrators are applied to an `IODE`.

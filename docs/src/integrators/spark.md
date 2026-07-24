```@meta
CurrentModule = GeometricIntegrators.SPARK
```

# Special Partitioned Additive Runge-Kutta Integrators

SPARK or *Special Partitioned Additive Runge-Kutta Integrators* are a family of integrators that have been introduced by [Laurent O. Jay](http://homepage.divms.uiowa.edu/~ljay/) for the integration of differential algebraic equations and in particular systems subject to holonomic and nonholonomic constraints [[Jay:1998](@cite), [Jay:2003](@cite), [Jay:2006](@cite), [Jay:2007a](@cite), [Jay:2007b](@cite)].
Recently, the idea of SPARK methods has been generalized and adapted to facilitate the integration of degenerate Lagrangian systems as well as Hamiltonian systems subject to Dirac constraints [[Kraus:2020](@cite)].

GeometricIntegrators.jl provides several flavours of such SPARK methods (*some are still experimental*):

| Integrator                          | Description                                                                                          |
|:------------------------------------|:-----------------------------------------------------------------------------------------------------|
| [`IntegratorHPARK`](@ref)           | Partitioned additive methods for Hamiltonian system subject to a general constraint $\phi(q,p) = 0$  |
| [`IntegratorVPARK`](@ref)           | Partitioned additive methods for Lagrangian system subject to a general constraint $\phi(q,p) = 0$   |
| [`IntegratorSPARK`](@ref)           | SPARK methods for general index-two differential algebraic equations                                 |
| [`IntegratorHSPARK`](@ref)          | Hamiltonian system subject to a general constraint $\phi(q,p) = 0$                                   |
| [`IntegratorHSPARKprimary`](@ref)   | Hamiltonian system subject primary constraint in the sense of Dirac                                  |
| [`IntegratorHSPARKsecondary`](@ref) | Hamiltonian system enforcing primary & secondary Dirac constraint                                    |
| [`IntegratorVSPARK`](@ref)          | Lagrangian system in implicit form subject to a general constraint $\phi(q,p) = 0$                   |
| [`IntegratorVSPARKprimary`](@ref)   | Degenerate Lagrangian system subject primary constraint in the sense of Dirac                        |
| [`IntegratorVSPARKsecondary`](@ref) | Degenerate Lagrangian system enforcing primary & secondary Dirac constraint                          |

These integrators are applied to either an `IDAE`, `HDAE` or `LDAE`.

## Construction

A SPARK method combines two Runge–Kutta tableaus: an *internal* pair
``(a_{ij}, b_i)`` on ``s`` stages that advances the dynamics, and a *projective*
(or *constraint*) pair ``(\tilde{a}_{ij}, \tilde{b}_i)`` on ``r`` stages that
enforces the algebraic constraint. Writing the vector field in the additive split
``f = f^{1} + f^{2} + f^{3}`` (Hamiltonian part, fibre-derivative part, constraint
part), the internal stages read

```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum_{j=1}^{s} a_{ij} V_{n,j} + h \sum_{j=1}^{r} \alpha_{ij} U_{n,j} , \\
P_{n,i} &= p_{n} + h \sum_{j=1}^{s} a_{ij} F_{n,j} + h \sum_{j=1}^{r} \alpha_{ij} G_{n,j} ,
\end{aligned}
```

the projective stages carry the Lagrange multipliers ``\Lambda_{n,i}`` and the
projective forces ``U_{n,i}``, ``G_{n,i}``, and the constraint is averaged over the
projective stages,

```math
0 = \sum_{j=1}^{r} \omega_{ij} \, \tilde{\Phi}_{n,j} + \omega_{i,r+1} \, \phi(q_{n+1}, p_{n+1}) ,
\qquad i = 1, \dots, r ,
```

with the last row of ``\omega`` selecting the constraint ``\phi(q_{n+1},p_{n+1})``
at the solution. The building blocks live in `RungeKutta.jl` (Gauß/Lobatto
tableaus) and in `src/spark/` (`gauss_ω_matrix`, `lobatto_ω_matrix`,
`lobatto_gauss_coefficients`, `get_lobatto_nullvector`); the concrete method
factories `SPARKGLRK`, `SPARKLob{ABC,ABD}`, `TableauVSPARK…`, `TableauHSPARK…`,
`SLRKLobatto…` assemble the internal/projective pairs.

The *degenerate variational* (`VSPARK`, `SLRK`) and *Dirac-constraint* (`HSPARK`)
families additionally enforce the **secondary** constraint
``\psi = \dot{p} - \nabla\vartheta(q)\cdot v = 0`` (the time derivative of the
primary constraint), following [[Kraus:2020](@cite)].

## Order and stability caveats

Verification against the source manuscripts and empirical convergence measurement
(see `test/verification/spark_convergence_tests.jl` and the SPARK-submodule pass in
the repository's verification report) established that several SPARK methods do not
reach — or do not converge to — their nominal order. These are **inherent
properties of the methods**, not implementation defects:

* A SPARK method cannot be simultaneously symplectic and enforce the constraint
  ``\phi(q_{n+1},p_{n+1}) = 0`` at the solution. Methods that do so (e.g. the
  Lobatto `IIIA-IIIB` / `IIIB-IIIA` `SPARK`/`HPARK` pairs) exhibit reduced order or
  divergence on the degenerate Lotka–Volterra system.
* The projection symplecticity conditions require ``R(\infty) = 1``. Since the
  implementation uses ``R(\infty) = (-1)^{s+1}``, the Gauß-based `SPARKGLVPRK` and
  `HPARKGLRK` methods drop from order ``2s`` to order ``2`` at ``s = 2``.
* Tableau pairs that coincide give a degenerate (singular) stage system, most
  visibly at ``s = 2`` (`VSPARK(SPARK…(2))`).

Methods that are *not* subject to these obstructions — `SLRK`, `SPARKGLRK`,
`SPARKLob{ABC,ABD}`, the symplectic-projection `VPARK`/`VSPARK` variants,
`VSPARKsecondary`, and `HSPARK` on Gauß / Lobatto `ABC`/`ABD` — reach exactly their
documented order (Gauß ``2s``, Lobatto ``2s-2``). The `HSPARKsecondary` family is
**experimental** and currently raises a `SingularException`.

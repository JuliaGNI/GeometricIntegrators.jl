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
| [`IntegratorSLRK`](@ref)            | Specialised Lobatto RK methods for degenerate Lagrangians, primary constraint at every stage         |

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

## Specialised Lobatto Runge–Kutta methods (`SLRK`)

`SLRK` is the one family in this submodule that is *not* a SPARK method in Jay's
sense. It is aimed specifically at degenerate Lagrangians of the form

```math
L (q, v) = \vartheta (q) \cdot v - H (q) ,
```

for which the fibre derivative ``p = \vartheta(q)`` cannot be inverted for ``v``.
The equations of motion are then an index-two DAE with the primary constraint
``\phi(q,p) = p - \vartheta(q) = 0``, and they preserve the *noncanonical*
symplectic two-form built from ``\Omega = \nabla\vartheta^{T} - \nabla\vartheta``.

### What makes it different

Jay's SPARK methods **average** the primary constraints over the internal stages,

```math
0 = \sum_{j=1}^{s} \omega_{ij} \, \phi (Q_{n,j}, P_{n,j}) , \qquad i = 1, \dots, s-1 ,
```

which is what breaks symplecticity. `SLRK` instead imposes the primary constraint
**individually at every stage** and averages the *secondary* constraints
``\psi = \dot{p} - \dot{q}\cdot\nabla\vartheta(q)`` in their place, appending
``\phi(q_{n+1},p_{n+1}) = 0`` as the last column of ``\omega``. There is a single
set of ``s`` stages — the projective stages coincide with the internal ones
(``\sigma = s``), so the multipliers ``\Lambda_{n,i}`` live on the same nodes as
the velocities ``V_{n,i}``.

### The scheme

```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum_{j=1}^{s} a^1_{ij} V_{n,j} + h \sum_{j=1}^{s} a^2_{ij} \Lambda_{n,j} , \\
P_{n,i} &= p_{n} + h \sum_{j=1}^{s} \bar{a}^1_{ij} F_{n,j} + h \sum_{j=1}^{s} \bar{a}^2_{ij} G_{n,j} , \\
0 &= \phi (Q_{n,i}, P_{n,i}) - \frac{d_i}{\bar{b}^1_i} \, \mu , \\
0 &= \sum_{j=1}^{s} \omega_{ij} \, \Psi_{n,j} + \omega_{i,s+1} \, \phi (q_{n+1}, p_{n+1}) , \\
0 &= \sum_{i=1}^{s} d_i \, V_{n,i} ,
\end{aligned}
```

with the update

```math
q_{n+1} = q_{n} + h \sum_{i=1}^{s} \big( b^1_{i} V_{n,i} + b^2_{i} \Lambda_{n,i} \big) ,
\qquad
p_{n+1} = p_{n} + h \sum_{i=1}^{s} \big( \bar{b}^1_{i} F_{n,i} + \bar{b}^2_{i} G_{n,i} \big) ,
```

and

```math
F_{n,i} = \frac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , \qquad
G_{n,i} = \nabla \vartheta (Q_{n,i}) \cdot \Lambda_{n,i} , \qquad
\Psi_{n,i} = F_{n,i} - V_{n,i} \cdot \nabla \vartheta (Q_{n,i}) .
```

The four tableaus are the struct fields `q` ``\to (a^1, b^1)``, `q̃` ``\to (a^2, b^2)``,
`p` ``\to (\bar{a}^1, \bar{b}^1)``, `p̃` ``\to (\bar{a}^2, \bar{b}^2)``. Every shipped
constructor sets `q̃ = q` and `p̃ = p`, so a method is determined by a single
conjugate-symplectic Lobatto pair.

The ``\omega`` matrix from `lobatto_ω_matrix(s)` is ``s \times (s+1)``; the system
``\omega\,[\Psi_{n,1},\dots,\Psi_{n,s},\phi(q_{n+1},p_{n+1})]^{T} = 0`` is
equivalent to the ``s-1`` Lobatto-IIIA-averaged secondary constraints together with
``\phi(q_{n+1},p_{n+1}) = 0``.

The Lobatto stage system is rank deficient by one in the ``V``-direction; the
multiplier ``\mu`` relaxes the primary constraint along the null vector ``d``
(`get_lobatto_nullvector`) and the extra condition ``\sum_i d_i V_{n,i} = 0``
removes the deficiency — the same construction as in `VPRK`. ``\mu`` appears in the
primary-constraint equation and nowhere else; adding it to the momentum-stage equation as
well makes the stage Jacobian singular at ``h = 1``, which is what finding S10 of
`VERIFICATION_REPORT.md` fixed.

### Available methods

| Constructor | ``(a^1, \bar{a}^1)`` pair | Order |
|:---|:---|:---|
| `SLRKLobattoIIIAB(s)`  | Lobatto IIIA / IIIĀ | ``2s-2`` |
| `SLRKLobattoIIIBA(s)`  | Lobatto IIIB / IIIB̄ | ``2s-2`` |
| `SLRKLobattoIIICC̄(s)`  | Lobatto IIIC / IIIC̄ | ``2s-2`` |
| `SLRKLobattoIIIC̄C(s)`  | Lobatto IIIC̄ / IIIC | ``2s-2`` |
| `SLRKLobattoIIID(s)`   | Lobatto IIID / IIID̄ | ``2s-2`` |
| `SLRKLobattoIIIE(s)`   | Lobatto IIIE / IIIĒ | ``2s-2`` |

Each ``\bar{X}`` tableau is the conjugate-symplectic partner of ``X`` computed by
`RungeKutta.symplectic_conjugate_coefficients`, so
``\bar{b}_i a_{ij} + b_j \bar{a}_{ji} = \bar{b}_i b_j`` and ``b = \bar{b}`` hold by
construction (verified numerically to ``3\times10^{-17}`` for ``s = 2,3,4``).

### Choice of gauge

The Lagrangian ``L = \vartheta(q)\cdot v - H(q)`` is unchanged, as a dynamical system,
by a gauge transformation ``\vartheta \to \vartheta + \nabla F``. The *discrete method*
is not, and for two of the six families the choice matters:

| | ``\mathrm{rank}(a^1)`` | ``\mathrm{rank}(\bar{a}^1)`` | usable when a component of ``\vartheta`` vanishes |
|:---|:---|:---|:---|
| `SLRKLobattoIIIAB(s)` | ``s-1`` | ``s-1`` | **no** — `SingularException` |
| `SLRKLobattoIIIBA(s)` | ``s-1`` | ``s-1`` | **no** — singular or divergent |
| `SLRKLobattoIIICC̄(s)` | ``s`` | ``s-1`` | yes |
| `SLRKLobattoIIIC̄C(s)` | ``s-1`` | ``s`` | yes |
| `SLRKLobattoIIID(s)`  | ``s`` | ``s`` | yes |
| `SLRKLobattoIIIE(s)`  | ``s`` | ``s`` | yes |

Lobatto IIIA has a zero first row and IIIB a zero last column, and each conjugate
partner inherits the defect — the `IIIAB`/`IIIBA` pairs are exactly the two for which
``a^1`` *and* ``\bar{a}^1`` are both rank deficient. On a one-form of full rank the
missing direction is supplied by the projection force ``G = \nabla\vartheta^{T}\Lambda``;
if a component of ``\vartheta`` vanishes identically then ``\nabla\vartheta^{T}`` loses
rank, ``G`` can no longer reach the corresponding momentum component, and the stage
system becomes singular at *every* step size.

Concretely, `GeometricProblems.LotkaVolterra2d` and `LotkaVolterra2dSingular` differ
only by the exact form ``d(q_1 q_2)``. `SLRKLobattoIIIAB` and `SLRKLobattoIIIBA` run on
the first and fail on the second; the other four give **identical** trajectories on
both, to ``\le 7\times10^{-16}`` over ``T = 1``. Prefer `SLRKLobattoIIID` or
`SLRKLobattoIIIE`, which are full rank in both blocks and also have the smallest
symplecticity defect.

### Requirements on the problem

!!! warning
    `SLRK` requires the `LDAEProblem`'s `f` to be the **full** force
    ``\partial L/\partial q = -\nabla H + \nabla\vartheta\cdot v``, because ``F`` is
    used directly both in the ``p`` update and as ``\dot{p}`` in the secondary
    constraint. This is the *opposite* convention to `VSPARKsecondary`, which
    expects `f = -∇H` and reconstructs ``\partial L/\partial q`` as `f + g(V)`.
    Passing a problem built for `VSPARKsecondary` to `SLRK` silently drops the
    ``\nabla\vartheta\cdot v`` term.

    `GeometricProblems.LotkaVolterra2d` ships both variants: use
    `ldaeproblem_slrk` with `SLRK` and `ldaeproblem` with `VSPARKsecondary`.

```julia
using GeometricIntegrators, GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d

prob = ldaeproblem_slrk([1.0, 1.0]; timespan = (0.0, 100.0), timestep = 0.1)
sol  = integrate(prob, SLRKLobattoIIIAB(3))
```

### What `SLRK` does and does not preserve

* **Primary constraint: exactly.** ``\max_n |\phi(q_n,p_n)| \approx 2\times10^{-15}``
  over ``10^3`` steps on Lotka–Volterra, for every constructor and every ``s``.
  This is the property the family exists for, and it is the one variational
  integrators lack.
* **Order: as documented.** ``2s-2``, confirmed empirically for all six families at
  ``s = 2, 3`` in `test/verification/spark_convergence_tests.jl`.
* **Symplecticity: approximate, not exact.** The source manuscript states exact
  preservation of the noncanonical two-form, but its proof leaves an uncontrolled
  term in the multiplier block — constraint ``\sum_i d_i V_{n,i} = 0`` has no
  counterpart for ``\Lambda``, so
  ``h \sum_j b^2_j (d_j/\bar{b}^1_j)\, \mathrm{d}\mu \wedge \mathrm{d}\Lambda_{n,j}``
  survives. Measured on Lotka–Volterra, the first Poincaré invariant ``\oint p\,dq``
  drifts secularly, by ``O(h^{p+1})`` per step, e.g. ``2.6\times10^{-4}`` after one
  step and ``1.6\times10^{-1}`` after 100 steps at ``h = 0.1`` for
  `SLRKLobattoIIIAB(2)`. Because `SLRK` holds ``\phi`` to round-off, its canonical
  ``\oint p\,dq`` and noncanonical ``\oint\vartheta(q)\cdot dq`` agree to every digit,
  so either one measures its symplecticity. The control is a genuinely symplectic
  variational integrator measured on the **canonical** invariant — `VPRKGauss(2)` and
  `VPRKGauss(3)` hold it at ``\sim 10^{-14}`` out to ``3\times10^{3}`` steps. Their
  *noncanonical* value is not a control and is not small (``10^{-6}`` to ``10^{-2}``),
  because a variational integrator leaves ``\{\phi = 0\}``. `SLRKLobattoIIID` and
  `SLRKLobattoIIIE` have the smallest defect. See the "Fifth pass" section of the
  repository's `VERIFICATION_REPORT.md`.
* **The symplecticity defect is gauge invariant.** Measured on the two
  gauge-equivalent Lotka–Volterra one-forms, the `∮p·dq` drift, the
  ``\|J^{T}\Omega J - \Omega\|`` defect and the energy drift agree digit for digit for
  every family that is usable on both. The defect is a property of the discretisation,
  not of how the Lagrangian is written — but see "Choice of gauge" above for which
  families are usable at all.

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

Note that reaching the documented order is a statement about *accuracy* only. As
described above, `SLRK` attains its full order and preserves the primary constraint
exactly, but is only approximately symplectic — the first bullet above (symplectic
*and* constraint-at-the-solution is impossible) applies to it too, since `SLRK` does
enforce ``\phi(q_{n+1},p_{n+1}) = 0`` through the last column of ``\omega``.

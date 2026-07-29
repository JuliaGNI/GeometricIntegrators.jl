# Integrator Verification Report

Verification of all integrators in **GeometricIntegrators.jl** (excluding the
SPARK submodule) for numerical correctness and consistency with the docstrings,
the documentation, and the geometric-numerical-integration literature.

Two complementary passes were used:

* **Static** — the implemented stage/update equations (`components!`,
  `residual!`, `update!`, `integrate_step!`) and the method trait/order
  declarations were cross-checked against the source docstrings, the
  `docs/src/integrators/*.md` pages, and the literature.
* **Dynamic** — a new numerical-verification test suite (`test/verification/`)
  empirically measures the convergence order of every method and the
  energy behaviour of the geometric integrators. Convergence orders were
  measured on the **nonlinear** pendulum / Lotka–Volterra problems where the true
  (nonlinear) order matters, because the harmonic oscillator is linear and lets
  some methods *super-converge* beyond their nominal order.

Concrete Butcher-tableau coefficients (Gauß/Radau/Lobatto nodes, weights, `aᵢⱼ`)
are defined in the external **RungeKutta.jl** package, not in this repository, so
several findings below have their root cause there; the dynamic tests exercise
those coefficients indirectly (a wrong tableau shows up as a wrong empirical
order).

---

## Summary of findings

| # | Severity | Method(s) | Issue | Root cause | Action |
|:--|:---------|:----------|:------|:-----------|:-------|
| 1 | correctness (**resolved**) | `KraaijevangerSpijker` | `order()` and the docs claimed 2, but the method is order **1** | Order *attribute* was wrong; the tableau genuinely satisfies only the order-1 conditions | Fixed upstream in RungeKutta.jl **v0.5.22** (`order()` now returns 1); `rk.md` order table corrected 2→1; verification test now asserts order 1 |
| 2 | correctness | `LobattoIIIF`, `LobattoIIIG` and their partitioned/VPRK variants | Documented "order 2s"; achieve order **2s−2** (order 2 for s=2) | "order 2s" construction not realised in the coefficients (RungeKutta.jl) | Recorded as `@test_broken`; reported |
| 3 | correctness | `LobattoIIIB(s)` as a **standalone ODE** integrator | Documented 2s−2, achieves order **1** | Singular `A` ⇒ degenerate stage system when used standalone | Recorded as `@test_broken`; documented as partitioned/variational-only |
| 4 | correctness (**fixed**) | `McLachlan4` | Documented 4th-order composition achieved only order **2** | Composition-*ordering* bug — each stage applied `φ` before its adjoint `φ*`; the √19 coefficient values themselves are the correct McLachlan (1995) values | **Fixed** (swap the `a`/`b` vectors in the `SplittingCoefficientsNonSymmetric` constructor); now converges at order 4 |
| 5 | correctness | `VPRKLobattoIIIAIIIB(3)` | Documented 2s−2 = 4, achieves order **2** (the plain PRK `LobattoIIIAIIIB(3)` reaches 4) | VPRK-specific | Recorded as `@test_broken`; reported |
| 6 | doc consistency | `VPRK` docstring | Symplecticity condition written `bᵢāᵢⱼ + bⱼaⱼᵢ = bᵢbⱼ` (unbarred), inconsistent with IRK/IPRK/`rk.md`/`vprk.md` | — | **Fixed** to general barred form |
| 7 | doc bug | `rk.md` partitioned & implicit-equation blocks | Stage-time index typos (`c_j` in `i`-sums; `∑ᵢ` where `∑ⱼ` meant) | — | **Fixed** |
| 8 | doc bug | `hpi.md` | Full-width digit `０` in the action-integral bound | — | **Fixed** |
| 9 | doc gap | `dvi.md`, `hpg.md` | Empty stub pages | — | **Written** (see Documentation below) |
| 10 | observation | `DVRK(Gauss(s))` | Full order 2s in class; **order reduction to s** when `ϑ` is given in a gauge with no vanishing components | Not a feature of degenerate Lagrangians — the method was being applied outside its hypothesis class | Corrected; in-class and out-of-class regimes both tested |
| 11 | observation | `DVIA`/`DVIB`/`CMDVI`/`CTDVI` | Accurate short-time, but do **not** converge cleanly to the ODE reference over long times | Low-order degenerate methods; needs maintainer interpretation | Fixed-step regression retained; reported |
| 12 | observation | `CGVI(Gauss(s))` | Position converges at order **2s−2**; `CGVI(Gauss(1))` is non-convergent (order 0) | — | Documented; order asserted for s ≥ 2 |

Finding 4 (McLachlan4) was an integrator-implementation bug in this repository
and has been **fixed** (see below). Findings 1–3, 5 and 10–12 concerning the
RK/Lobatto coefficients ultimately originate in **RungeKutta.jl** or are inherent
properties of the methods; they are **not** integrator-implementation bugs here.
The integrator
machinery itself is correct: every collocation/standard method
(ERK, `Gauss`, `RadauIIA`, `LobattoIIIA/C/D/E`, `PartitionedGauss`,
`LobattoIIIAIIIB(2)`, Strang/TripleJump/SuzukiFractal, PMVI, DEL, HPI, projected
Gauß) reaches exactly its documented order.

---

## Verified-correct methods (dynamic order confirmation)

Empirical convergence orders match the documented order for:

* **Explicit RK**: `ExplicitEulerRK`(1), `ExplicitMidpoint`/`Heun2`/`Ralston2`/`Runge2`/`SSPRK2`(2), `Heun3`/`Kutta3`/`Ralston3`/`SSPRK3`(3), `RK4`/`RK438`(4).
* **DIRK**: `CrankNicolson`(2), `Crouzeix`(3), `QinZhang`(2).
* **Fully implicit RK**: `ImplicitEulerRK`(1), `ImplicitMidpoint`(2), `SRK3`(4), `Gauss(s)`(2s), `RadauIA`/`RadauIIA`(2s−1), `LobattoIIIA`/`IIIC`/`IIID`/`IIIE`(2s−2). `RadauIB`/`RadauIIB` confirmed order 2s−1 (their order-4 reading on the harmonic oscillator is linear super-convergence).
* **Partitioned RK**: `SymplecticEulerA`/`B`(1), `PartitionedGauss(s)`(2s), `LobattoIIIAIIIB`/`LobattoIIIBIIIA`(2s−2).
* **Splitting/composition**: `LieA`/`LieB`(1), `Strang`/`McLachlan2`(2), `McLachlan4`/`TripleJump`/`SuzukiFractal`(4, `McLachlan4` after the fix below).
* **Variational**: `PMVImidpoint`/`PMVItrapezoidal`(2), `DiscreteEulerLagrange` (Midpoint/Trapezoidal)(2), `VPRKGauss(s)`(2s), `VPRKLobattoIIIAIIIB(2)`(2).
* **Galerkin**: `CGVI(Gauss(s))` at order 2s−2 for s ≥ 2.
* **Degenerate variational**: `DVRK(Gauss(s))` at order 2s on every problem given in a gauge satisfying its hypothesis (`ϑ_μ = 0` for `μ > d/2`) — the degenerate pendulum, `LotkaVolterra2dSingular`, `MasslessChargedParticleSingular`; order s outside that class — see finding 10.
* **Hamilton–Pontryagin**: `HPImidpoint`/`HPItrapezoidal`(2).
* **Projected**: `PostProjection`/`MidpointProjection`/`SymmetricProjection` with `Gauss(s)`(2s), each conserving the Hamiltonian to machine precision.

---

## Detail on the order deficiencies

**KraaijevangerSpijker (finding 1).** Its tableau (from RungeKutta.jl) has
`c = [0.5, 1.5]`, `b = [-0.5, 1.5]`, giving `Σᵢbᵢ = 1` but `Σᵢbᵢcᵢ = 2 ≠ 1/2`.
The order-2 condition is therefore violated at the coefficient level, so the
method is genuinely order 1, whereas `order(KraaijevangerSpijker())` and the
`rk.md` table previously stated order 2. **Resolved:** RungeKutta.jl v0.5.22
corrects the order attribute to 1 (the coefficients are unchanged and remain a
valid order-1 method); the `rk.md` order table has been updated to 1 and the
verification test now asserts order 1.

**LobattoIIIF / LobattoIIIG (finding 2).** Documented as order `2s` (the Fangzong
construction, [Fangzong:2016]). Measured order is `2s−2` in every context tested
— standalone (`LobattoIIIF(2)`→2), as symplectic partitioned pairs
(`LobattoIIIFIIIF̄(2)`, `LobattoIIIGIIIḠ(2)`→2), and as VPRK methods
(`VPRKLobattoIIIF̄(2)`, `VPRKLobattoIIIG(2)`→2; `LobattoIIIG(3)`→4). The order-2s
property is not realised by the coefficients.

**LobattoIIIB standalone (finding 3).** `LobattoIIIB(2)` has a singular
coefficient matrix (`A = [[½,0],[½,0]]`), so applied directly to an ODE its stage
system is degenerate and it drops to order 1. This method is intended for use in
symplectic *partitioned* pairs (`LobattoIIIAIIIB`, order 2s−2, verified correct)
and as a VPRK component, not as a standalone ODE integrator.

**McLachlan4 (finding 4) — verified against the literature and fixed.** The
√19 coefficients in the code exactly match the method's own docstring, and they
are the correct McLachlan (1995) 4th-order "minimum error" coefficients. The
order-2 behaviour was **not** a coefficient error but a *composition-ordering*
bug. This was confirmed with an independent matrix (BCH) model on non-commuting
generators `A`, `B`, comparing one step of the composition against
`exp(h(A+B))` (local-error slope = order + 1):

* Composing a **first-order method `φ` and its adjoint `φ*`** with each stage as
  `φ` *then* `φ*` (what `SplittingCoefficientsNonSymmetric` produces from
  `(a, b) = (α, β)`) → slope 3 ⇒ **order 2**.
* The *same* coefficients with each stage as `φ*` *then* `φ` (i.e. swapping `a`
  and `b`), or composed on a second-order **symmetric Strang** base → slope 5 ⇒
  **order 4**.

Corroborating evidence for the original deficiency: the existing accuracy test
(`test/integrators/splitting_integrators_tests.jl`) used `< 5E-4` for
`McLachlan4`, looser than the `< 1E-4` used for the *second*-order `McLachlan2`.

**Literature cross-check (McLachlan 1995, Table 2).** The code's `McLachlan4`
vectors are `a = [c₁,…,c₅]` and `b = reverse(a) = [d₁,…,d₅]`. Starting from the
paper's order-4 *type-S, m=5* minimum-error parameters
(`b₁=2/5, b₂=−1/10, a₁=(14−√19)/108, a₂=(20−7√19)/108`) and solving the
determining equations (2.6) `d₁=a₁`, `dᵢ+cᵢ=bᵢ`, `dᵢ+cᵢ₋₁=aᵢ` for the actual
composition coefficients `cᵢ,dᵢ` reproduces the code's vectors **exactly** (to
machine precision). So the coefficient values are confirmed correct against the
primary source; only the composition ordering was wrong. (The fix produces the
`A↔B`-relabelled composition, which for a symmetric split is an equivalent
4th-order method with the same error constant.)

**Fix.** In `coefficients(::McLachlan4, …)`, the constructor call was changed from
`SplittingCoefficientsNonSymmetric(:McLachlanSplitting, 4, a, a[end:-1:1])` to
`SplittingCoefficientsNonSymmetric(:McLachlanSplitting, 4, a[end:-1:1], a)`,
which composes each stage's adjoint flow first. The real integrator now converges
at order 4.02 on the harmonic oscillator, and the existing splitting tests still
pass. (`McLachlan2` uses the same coefficient pattern but is order-correct at 2
regardless of the stage ordering, so it was left unchanged.)

**VPRKLobattoIIIAIIIB(3) (finding 5).** Reaches order 2 rather than 2s−2 = 4,
although the plain partitioned `LobattoIIIAIIIB(3)` reaches 4. The reduction is
specific to the VPRK implementation of this pair at s = 3.

## Detail on the observations

**DVRK order and the gauge of ϑ (finding 10, corrected).** The earlier reading of
this finding — "order reduction to s is a known feature of symplectic RK methods
on degenerate Lagrangians" — was wrong. The reduction is not a property of
degenerate Lagrangians; it is what happens when `DVRK` is applied *outside its
hypothesis class*.

`DVRK` requires `L = ϑ(q)⋅q̇ − H(q)` with `d/2` components of `ϑ` vanishing
identically. `ϑ` is only defined up to an exact one-form, and the gauge decides
whether a given problem satisfies that hypothesis. Every DVRK test previously in
the suite used `GeometricProblems.LotkaVolterra2d`, whose gauge is
`ϑ = (q₂ + log q₂/q₁, q₁)` with `ϑ₂ ≠ 0` — outside the class. Measured orders
(T = 1, reference `Gauss(8)` on the ODE):

| Problem | Gauge | s = 1 | s = 2 | s = 3 |
|---|---|---|---|---|
| `LotkaVolterra2dSingular`         | `ϑ₂ = 0` (in class) | 2.00 | 4.00 | 5.98 |
| `LotkaVolterra2d`                 | `ϑ₂ = q₁ ≠ 0`       | 1.01 | 2.00 | 3.01 |
| `LotkaVolterra2dSymmetric`        | both ≠ 0            | 0.99 | 2.00 | 3.01 |
| `MasslessChargedParticleSingular` | `A₂ = 0` (in class) | 2.00 | 4.00 | 5.98 |
| `MasslessChargedParticle`         | both ≠ 0            | 1.00 | 2.00 | 3.00 |

Within each block the rows describe the *same dynamics*, differing only by a gauge
transformation of `ϑ`. In class, `DVRK(Gauss(s))` attains the full order 2s; out of
class it remains convergent at order s. Both regimes are now covered by
`test/verification/dvi_convergence_tests.jl`, and `docs/src/integrators/dvi.md`
states the hypothesis and the gauge dependence explicitly.

The split is between gauges, not between regular and degenerate Lagrangians: the
pendulum `IODEProblem`, previously recorded here as a *regular* Lagrangian on which
`DVRK` reaches 2s, is in fact the degenerate phase-space form with
`ϑ = (m l² q₂, 0)` on `q = (x, p)` — that is, in class. Every 2s measurement in the
suite is in class and every order-s measurement is out of class.

The invertibility hypothesis on `∂ϑ_μ/∂q^ν` (`μ ≤ d/2 < ν`) is a separate and
weaker condition, and does not detect the out-of-class cases: for
`LotkaVolterra2d`, `∂ϑ₁/∂q² = 1 + 1/q₁q₂ ≠ 0` is invertible while `ϑ₂ ≠ 0`.

**Low-order degenerate integrators (finding 11).** `DVIA`, `DVIB`, `CMDVI`,
`CTDVI` reproduce the established short-time accuracy on Lotka–Volterra
(`T = 0.1`), but over longer times (`T = 1`) their error against the ODE
reference does not decrease under step refinement. This may be an inherent
limitation of these first/second-order degenerate integrators on this problem
rather than a defect; it is flagged for maintainer interpretation. The tests
retain the established short-time accuracy checks.

**CGVI order (finding 12).** The continuous Galerkin integrator converges at
order `2s−2` in the position (`Gauss(2)`→2, `Gauss(3)`→4, `Gauss(4)`→6);
`CGVI(Gauss(1))` does not converge (order 0). No order was previously claimed in
the docs, so this is documentation rather than a discrepancy.

---

## Fixes applied

* `src/integrators/vi/vprk_methods.jl` — `VPRK` docstring symplecticity condition
  corrected from the unbarred `bᵢāᵢⱼ + bⱼaⱼᵢ = bᵢbⱼ` to the general barred form
  `bᵢāᵢⱼ + b̄ⱼaⱼᵢ = bᵢb̄ⱼ`, matching `integrators_irk.jl`, `integrators_iprk.jl`,
  `docs/src/integrators/rk.md`, and `docs/src/integrators/vprk.md`.
* `docs/src/integrators/rk.md` — fixed stage-time and summation-index typos in the
  partitioned-equation and implicit-equation blocks (`c_j`→`c_i` in the `i`-sums;
  `∑ᵢ`→`∑ⱼ` in the `Pₙ,ᵢ` stage equation).
* `docs/src/integrators/hpi.md` — replaced the full-width digit `０` with `0` in
  the action-integral lower bound.
* `src/integrators/splitting/splitting_methods.jl` — fixed the `McLachlan4`
  composition ordering (finding 4): the `SplittingCoefficientsNonSymmetric`
  constructor now receives the coefficient vectors swapped so that each stage's
  adjoint flow is composed first, restoring the documented 4th order. The
  coefficient values are unchanged (verified correct against the literature).

## Non-findings (checked, correct)

* `ExplicitEuler`/`ImplicitEuler` aliases (referenced by `rk.md`) **are** defined
  and exported, so the documentation cross-references resolve.
* Stage-time convention `sol.t + h·(cᵢ − 1)` in `components!` is correct:
  `sol.t` is `t_{n+1}`, so the argument equals `t_n + h·cᵢ`, consistent with
  `initial_guess!` (`history[1].t + h·cᵢ`).
* `RadauIB`/`RadauIIB` are order 2s−1 as documented.

---

## Documentation added

* `docs/src/integrators/dvi.md` — full Degenerate Variational Integrators page
  (degenerate Lagrangian `L = ϑ(q)·q̇ − H(q)`; the `DVIA`/`DVIB`/`CMDVI`/`CTDVI`
  Euler/midpoint/trapezoidal methods; the `DVRK` scheme and symplecticity
  theorem), following the CGVI/VPRK style. References [Ellison:2018],
  [Kraus:2019], [Kraus:2017].
* `docs/src/integrators/hpg.md` — Hamilton–Pontryagin–Galerkin framework page
  (continuous HP principle, broken-Galerkin discretisation, Type I–IV continuity
  constraints / numerical fluxes). Marked clearly as a not-yet-implemented
  framework; re-enabled in `docs/make.jl`. References [Yoshimura:2006a],
  [Yoshimura:2006b], [BouRabee:2009], [Leok:2011].
* `docs/src/GeometricIntegrators.bib` — added `Yoshimura:2006a`, `Yoshimura:2006b`,
  `BouRabee:2009`, `Tyranowski:2014linear`, `Ellison:2018`, `Kraus:2019`.

## Tests added

New suite under `test/verification/`, wired into `test/runtests.jl`:

* `verification_utilities.jl` — `estimate_convergence_order`,
  `test_convergence_order`, `test_energy_behaviour`, `loglog_slope`.
* Per-family files: `rk`, `prk`, `splitting`, `variational`, `galerkin`, `dvi`,
  `hpi`, `projection` convergence tests. Confirmed order deficiencies (findings
  1–5) are recorded with `@test_broken` at the *documented* order, so that a
  future coefficient fix is detected automatically.

## Methods added

The package previously provided splitting/composition methods only up to 4th
order. Two higher-order symmetric composition methods from the literature
(Yoshida 1990; tabulated in McLachlan 1995, Table 2) were added in
`src/integrators/splitting/splitting_methods.jl` and registered in
`src/integrators/method_list.jl` (which auto-exports them and lists them in the
generated methods table):

* **`Yoshida6`** — 6th-order symmetric composition (SS, m = 7, Yoshida "solution A").
* **`Yoshida8`** — 8th-order symmetric composition (SS, m = 15, Yoshida "solution D").

Both use the same `SplittingCoefficientsSS` framework as `TripleJump`/`SuzukiFractal`.
Their coefficients were verified to high precision: a BigFloat matrix-exponential
order check (composing a symmetric Strang base) gives clean order 6 and 8 down to
`Δt = 0.025`, and the real integrator reaches order 5.95 / 8.04 on the harmonic
oscillator. (During transcription a ~1e-11 error in the 6th-order `w₁` was caught
by this high-precision cross-check — the low-precision test could not distinguish
it — and corrected against the value already recorded in `docs/src/integrators/splitting.md`.)
Convergence tests for both were added to `test/verification/splitting_convergence_tests.jl`
and type/order checks to `test/methods/splitting_methods_tests.jl`.

## Dead / disabled code (noted, not verified)

Present on disk but not compiled into the `Integrators` module, and therefore
outside the verification scope:
`src/integrators/dgvi/*` (DGVI, commented out), `src/integrators/vprk/*` (an
alternative VPRK implementation, not included), `src/integrators/VPRK.jl` (legacy
module), and the commented-out `rk/integrators_midpoint_implicit.jl`,
`integrators_srk_implicit.jl`, `integrators_flrk.jl`, `pglrk_integrators.jl`.
The CGVI documentation page remains disabled in `docs/make.jl` despite CGVI being
implemented; re-enabling it is left to the maintainers.

---

# Second pass — test-tolerance audit & disabled-test re-enablement

A follow-up pass over the **test suite** with two goals: (1) tighten error
tolerances that were more than an order of magnitude looser than the error the
test actually achieves, and (2) verify every disabled test of an existing
integrator — re-enable it where the integrator works, set an honest tolerance
from the measured error, and record the genuinely-broken cases as `@test_broken`
with a root cause. SPARK was included in this pass. Every tolerance below was set
empirically (measure the error, then take the smallest round `{1,2,4,5,8}×10ᵏ`
bound ≥ ~2× the measurement).

## Tolerances tightened (measured error ≫ 10× below the old bound)

* `test/integrators/rk_integrators_tests.jl` — `Gauss(3)` on the harmonic
  oscillator (ODE and partitioned): `1E-11 → 2E-12` (measured 7.5E-13).
* `test/integrators/galerkin_integrators_tests.jl` — `CGVI`: `1E-7 → 8E-13`
  (measured 3.7E-13).
* `test/integrators/ensemble_integrators_tests.jl` — ODE ensemble `Gauss(2)`:
  `8E-6 → 8E-8` (measured 3.5E-8).
* `test/projections/projections_implicit_tests.jl` — the Post/Midpoint/Symmetric
  projection `Gauss(1..4)` ladders (Lotka–Volterra): e.g. `4E-4 → 2E-6`,
  `8E-7 → 2E-11`, `4E-11 → 2E-15`, `4E-14 → 2E-15`.
* `test/projections/projections_vprk_tests.jl` — restored the previously-relaxed
  tighter values recorded in the inline `#…` comments (measurements confirm they
  hold): PostProjection `VPRKGauss(1) 4E-6 → 1E-6`, `(2) 8E-7 → 1E-11`,
  `(3) 4E-11 → 2E-15`; SymplecticProjection `VPRKGauss(2/3)` likewise.

The remaining direct-error families (explicit/implicit/DIRK RK, splitting, DVI,
HPI, PMVI/VPRK, SLRK, VPARK, most VSPARK/HSPARK) were already tuned to ~2–5× the
measured error and were left unchanged.

## Disabled tests re-enabled

* **Momentum (`p`) error checks** (RK partitioned, VPRK) — were disabled because
  `relative_maximum_error(sol, ref).p` evaluated `0/0 = NaN`: on the harmonic
  oscillator the reference momentum has `p(0) = 0`, and the metric normalised per
  timestep by the momentum amplitude. The momentum is in fact accurate. This was a
  bug in **GeometricSolutions**, fixed upstream in **v0.6.5**; the package's
  compat and the test manifest were updated to that release, and the `p`-checks
  were re-enabled calling `relative_maximum_error(...).p` (and
  `relative_maximum_error(sol.p, pref.p)`) directly, with tolerances set from
  measurement. (An interim `gerr` global-amplitude-normalised helper used before
  the upstream fix has been removed.)
* **Exact-equality "TODO: Reactivate!" checks** — in `rk_integrators_tests.jl`
  the ODE↔PODE↔HODE conversion results are bit-identical, so the `==` checks were
  restored (replacing the looser `≈ atol=1E-15`). In
  `splitting_integrators_tests.jl` the direct-vs-`Composition` results differ at
  ~1E-16 (floating-point non-associativity), so `==` does **not** hold; the
  approximate `< 1E-15`/`2E-15` checks were kept and the stale TODOs replaced with
  an explanatory note.
* **Projections** — `PostProjection(RK4())` re-enabled in `projections_tests.jl`
  (`< 4E-7`). `VPRKpVariationalQ`/`VPRKpVariationalP` re-enabled in
  `projections_vprk_tests.jl` with the stale `reference_solution` replaced by
  `ref.q`. `VPRKpSecondary` `show` re-enabled in `test_show.jl` (construction only).

## Disabled SPARK tests re-enabled / analysed (`test/spark/spark_integrators_tests.jl`)

About 30 previously-commented SPARK/VSPARK/HSPARK cases were re-enabled with honest
tolerances after confirming they run and converge, e.g.
`SPARKLobattoIIIAIIIB(3/4)`, `VSPARK(SPARKLob{ABC,ABD,attoIIIAIIIB,attoIIIBIIIA})`
for `s = 3,4`, and `TableauVSPARKLobattoIII{AIIIB,BIIIA}pSymmetric`. The
`SPARK/HSPARK(GLRKLobatto…)` families were re-enabled but exhibit **order
reduction** (accuracy plateaus at ~2E-4 regardless of `s`), which is documented in
the test comments.

21 cases were recorded as `@test_broken` with a root cause:
* **Singular stage system** at `s = 2` — `VSPARK(SPARKGLRK/SPARKLob…(2))` and
  `HSPARK(SPARKLobattoIII{AIIIB,BIIIA}(2/3/4))` raise `SingularException`.
* **Divergence / excessive solver iterations** — `TableauHPARKLobattoIII{AIIIB,BIIIA}(2/3)`
  (errors 0.15–20), `SPARKLobattoIIIBIIIA(2/3/4)`, `SPARKLobattoIIIAIIIB(2)`.
* **Order reduction** — `TableauVSPARKLobattoIIIBIIIApSymmetric(3)` (no better than `s = 2`).

The **HSPARK-secondary** testset (`TableauHSPARKLobatto…` on `hdae`) raised
`FieldError: type NamedTuple has no field q̇`. Root cause: `initial_guess!` in
`src/spark/integrators_hspark_secondary.jl` built the extrapolation NamedTuple
with fields `v`/`f`, whereas `solutionstep!` (the Hermite extrapolation) writes
the solution derivatives `sol.q̇`/`sol.ṗ` — every other SPARK integrator already
uses `q̇=`/`ṗ=` (sourced from the velocity/force caches `V`/`F`). **Fixed** by
renaming those two NamedTuple fields to `q̇`/`ṗ` (values unchanged), which
resolves the interface incompatibility. The testset was re-enabled; the methods
now run through the solver but still fail for two separate, pre-existing SPARK
numerical issues, recorded as `@test_broken`:
`TableauHSPARKLobattoIII{AB,BA,D,E}` raise a `SingularException` (singular stage
system at all orders), and `TableauHSPARKGLRKLobattoIII{AB,BA,D,E}` raise a
`BoundsError` (an `s×s` coefficient matrix is indexed at `[1, s+1]` inside
`components!`). Both are independent of the field-name fix.

## Modernisation attempts — blockers

The following disabled tests target code that was removed or predates the
method-based rearchitecture; reviving them is a re-implementation effort, not a
test port. They are left disabled with explanatory comments:

* **FLRK** (Formal Lagrangian RK) — the `FLRK` method is commented out in
  `method_list.jl` and the integrator source (`integrators_flrk.jl`) was removed.
* **PGLRK** — the integrator (`IntegratorPGLRK`, `pglrk_integrators.jl`) was
  removed, and `CoefficientsPGLRK(s)` is bit-rotted against current dependency
  APIs (unqualified `Legendre`; `QuadratureRules` `nodes`/`weights` now return
  matrices; an inner-constructor field-type mismatch remains even after those two
  fixes). The `IntegratorVPRKpTableau` path lives only in the uncompiled legacy
  `vprk/` module.

  > **Superseded by the fourth pass below.** All three families have since been
  > modernised. Two of the claims above are wrong and are corrected there: the
  > removed integrator sources are recoverable from `38c64da0^` (and the `FLRK`
  > *method* struct from `2d927d9f^:src/methods/flrk.jl`), the commented include
  > named a file (`pglrk_integrators.jl`) that never existed, and
  > `QuadratureRules` `nodes`/`weights` return `Vector`s, not matrices — the
  > 0-based `leg_basis[nodes, j-1]` indexing was always correct.
* **DGVI** — `src/integrators/dgvi/*` use the superseded `Parameters` /
  `update_params!` architecture; a full port to `GeometricIntegrator` is required.

  > **Superseded by the fourth pass below.** All five variants have since been
  > ported. The port was smaller than this estimate suggests: the numerics had not
  > regressed since `6245e4f3`, only the surrounding architecture.
* **Old-API RK tests** — the `AbstractIntegrator(ode, Tableau…())` /
  `IntegratorIRK(…; exact_jacobian=true)` / block-Jacobian tests use constructors
  removed in the rearchitecture. The integrator-type checks are already covered by
  `test/methods/runge_kutta_methods_tests.jl`; the block-Jacobian test has no
  current equivalent.

## Additional findings (not fixed)

* **`order(VPRKGauss(s))` / `order(VPSRK3())` metadata is wrong.** `VPRK(Gauss(1))`
  and `VPRK(VPRKGauss(1))` have **identical** Butcher coefficients, but the
  VPRK-alias tableau reports the wrong order: `order(VPRKGauss(1)) = 1` (vs
  `order(Gauss(1)) = 2`) and `order(VPSRK3()) = 3` (vs `order(SRK3()) = 4`). The
  root cause is RungeKutta.jl's `SymplecticPartitionedTableau` (the `o` field of
  the symplectic-conjugate partitioned tableau), analogous to findings 1–3/5. The
  `test/verification` convergence tests pass **explicit** expected orders, so they
  are unaffected; the disabled `VPRK(X) == VPRK(Y)` equivalence checks in
  `test/methods/vprk_methods_tests.jl` legitimately fail on this (and on the
  tableau `name`) metadata and are left disabled with a note.
* **Several VPRK projection wrappers are non-functional.** `VPRKpInternal`,
  `VPRKpSecondary`, `VPRKpVariational` (and `VPRKpLegendre`) are exported but
  `integrate` throws `MethodError: no method matching initial_guess!(…)` for their
  solution-step layout — the projection integrators are not fully wired into the
  current architecture. Only `VPRKpVariationalQ`/`VPRKpVariationalP` and the
  standard/symplectic/midpoint/symmetric projections work. Left disabled and flagged.
* **`MidpointProjection(RK4())` / `SymmetricProjection(RK4())`** error in
  `initial_guess!` (explicit RK4 supplies no vector field for the projection
  guess); left disabled.

---

# Third pass — spurious solver-iteration warnings in variational / Hamilton–Pontryagin convergence tests

The `variational` and `hpi` convergence suites printed
`Warning: Solver took 1000 iterations.` at the finest timestep of `steps(10, 4)`
(Δt = 1/160). Only the **position–momentum** (`PMVImidpoint`, `PMVItrapezoidal`)
and **Hamilton–Pontryagin** (`HPImidpoint`, `HPItrapezoidal`) integrators are
affected; DEL, all VPRK variants and CGVI are clean. On the finest-Δt run a
handful of the 160 per-step solves hit the cap: `PMVImidpoint` 2,
`PMVItrapezoidal` 3, `HPImidpoint` 3, `HPItrapezoidal` 2.

**Not a solver-quality problem.** The nonlinear solver is SimpleSolvers.jl's
Newton (full Newton, `refactorize = 1`) with a `Backtracking` line search. The
warning (`SimpleSolvers/…/nonlinear_solver_status.jl`) fires when
`iterations ≥ warn_iterations`, and since `warn_iterations = max_iterations = 1000`
it means the solve hit the cap without converging. On the affected cases the
warning count is **identical** for Newton+Backtracking (default), Newton+Static,
Newton+StrongWolfe and **DogLeg**; Bisection is worse and Quadratic much worse.
Raising `max_iterations` to 10 000 does not help — the solves are genuinely
**stalled**, not slow. In every configuration the empirical order stays exactly 2.

**Root cause.** The framework default `f_abstol = 8eps() ≈ 1.78e-15`
(`GeometricIntegratorsBase.default_options`) is essentially machine precision. At
a few phase points the residual stagnates just above ~1e-15 and can never meet
it, so Newton spins to the cap. This is a tolerance-stagnation artefact, not a
defect in the integrators or the solver.

**Fix (test-only).** The `variational` and `hpi` convergence suites now pass
`integrate_options = (min_iterations = 1, f_abstol = 4e-15)` to the four affected
integrator calls. A sweep locates the threshold precisely at `4e-15`: `3e-15`
still warns (once each for `PMVImidpoint` / `HPImidpoint`), while `4e-15` removes
all warnings for all four methods with order = 2 and all 5 points valid — still
far tighter than the `f_reltol ≈ 1.5e-8` default and only ~2.3× the
machine-precision default. `min_iterations = 1` is repeated because passing any
solver option replaces the whole `default_options` NamedTuple. The
`test_convergence_order` / `estimate_convergence_order` helpers in
`verification_utilities.jl` gained an `integrate_options` keyword (default empty,
backward-compatible) to thread these options through to `integrate`.

## SPARK suite (`test/spark/spark_integrators_tests.jl`)

The SPARK suite emitted ~35 warnings per run, of two kinds — `Solver took 1000
iterations.` (nonlinear solver hit its cap) and `Backtracking line search did not
satisfy the sufficient decrease condition within 1000 iterations.` Unlike the
PMVI/HP case these are **mostly not** a uniform benign artefact; they split three
ways:

* **Genuinely divergent / order-broken methods, already recorded as `@test_broken`**
  — the bulk of the noise. `SPARKLobattoIIIBIIIA(2/3/4)` (errors 2.96 / 0.58 /
  7.3e-3, up to 9 solver + hundreds of line-search warnings), `SPARKLobattoIIIAIIIB(2)`
  (0.107), `TableauVSPARKLobattoIIIBIIIApSymmetric(3)` (order-reduced), and the
  `TableauHPARKLobattoIII{AIIIB,BIIIA}(2/3)` cases (errors 0.15–20). The warnings
  are a *correct symptom* of divergence; a different solver (DogLeg), line search,
  or a higher iteration cap does not help (DogLeg makes matters worse or raises a
  `NonlinearSolverException`).
* **Singular / out-of-bounds `@test_broken` methods** — abort with
  `SingularException` / `BoundsError`; most warn zero times, except
  `VSPARK(SPARKLobABD(2))`, which emits a few warnings before the singular solve.
* **Benign warners among the passing cases** — `VSPARK(SPARKLobABD(4))` (1 solver +
  204 line-search warnings) and `VSPARK(SPARKLobABC(3))` (19 line-search warnings,
  solver converges). Both pass with full accuracy.

**Root cause & fix.** SPARK overrides `default_options` with an even tighter set
(`min_iterations=1, x_suctol=2eps(), f_abstol=8eps(), f_suctol=2eps()`,
`src/spark/abstract.jl`). For `VSPARK(SPARKLobABD(4))` the same machine-precision
`f_abstol` stall as PMVI/HP applies; relaxing it to `4e-15` (keeping the other
SPARK defaults, threaded via the module-level `SPARK_RELAXED` named tuple) makes
the solve converge and removes **both** warning kinds, error unchanged (5.2e-12).

For every other noisy case, tolerance/solver tuning does **not** cleanly help —
the divergent methods genuinely fail, and the line-search warning cannot be
silenced through the public API (the line search builds its own `Options` and
never receives a `verbosity` kwarg passed through `integrate`; `warn_iterations=0`
suppresses only the solver-cap warning). These calls are therefore wrapped in a
`muffle(f) = with_logger(f, NullLogger())` helper that suppresses log output for
that one integration. `muffle` changes only logging, so the measured errors and
the `@test_broken` status are unaffected. The suite now runs warning-free with the
same 133 pass / 45 broken tally.

---

# Fourth pass — SPARK submodule (`src/spark/`)

The first three passes covered the whole package **except** the SPARK submodule.
This pass completes the verification: a static + dynamic check of every SPARK
integrator family (`SPARK`, `VPARK`, `HPARK`, `VSPARK`/`primary`/`secondary`,
`HSPARK`/`primary`/`secondary`, `SLRK`) and a root-cause diagnosis of why several
SPARK methods do not converge.

## Method and reference basis

SPARK methods are largely **original** and not documented in the standard
geometric-numerical-integration literature. Verification used the two draft
manuscripts *SPARK Methods for Degenerate Lagrangian Systems* and *SPARK Methods
for Hamiltonian Systems Subject to Dirac Constraints*. Those manuscripts are
**unfinished**: they prove the *symplecticity* conditions but defer all
order/convergence analysis to Jay's cited papers ("order and convergence will be
discussed elsewhere"). Concrete order numbers therefore come from the code's
tableau `o` fields (Gauß `2s`, Lobatto-pair `2s−2`, SLRK-Lobatto `2s−2`) and from
empirical measurement. The manuscripts do, however, make three qualitative
predictions that explain almost all of the observed non-convergence:

1. **Symplectic ⇒ constraint-at-solution is impossible.** (Degenerate paper,
   Sec. 3.4 remark.) A SPARK method cannot be symplectic *and* enforce
   `φ(qₙ₊₁,pₙ₊₁)=0`; such methods "show reduced order of convergence or even
   divergence."
2. **`R(∞) ≠ 1` breaks the projection symplecticity conditions.** The code passes
   `R∞ = (-1)^(s+1)` explicitly, so odd/even stage counts flip whether the
   conditions hold.
3. **Coinciding tableau pairs are singular.** The velocity `V` and the multiplier
   `Λ` are indistinguishable when the two coefficient blocks coincide, giving a
   degenerate (singular) stage system — most visibly at `s = 2`.

## Dynamic verification added

A new convergence suite `test/verification/spark_convergence_tests.jl` (wired into
`test/runtests.jl`) measures the empirical order of every family on the degenerate
Lotka–Volterra system in its `IDAE`/`PDAE`/`LDAE`/`HDAE` formulations, referenced
against `Gauss(8)` on the ODE form, at `T = 1` (longer than the fixed-step
`spark_integrators_tests.jl`, `T = 0.1`, so the true asymptotic behaviour shows).
It asserts the documented order for the working methods (64 checks) and records the
confirmed deficiencies as `@test_broken` at the documented order (11). Empirically
confirmed orders:

* **Full documented order:** `SLRKLobattoIII{AB,BA,D,E}` (2s−2); `SPARKGLRK(s)`
  (2s); `SPARKGLVPRK(1)` (2); `SPARKLob{ABC,ABD}(s)` (2s−2); `TableauGausspSymplectic(2)`,
  `TableauLobattoIII{AIIIB,BIIIA}pSymplectic(3)`, `TableauVSPARKGLRKp{Midpoint,Symplectic,Symmetric}(2)` (4);
  `VSPARK(SPARKLobattoIIIAIIIB(3))` (4); `TableauVSPARKLobattoIIIAB(s)` and
  `TableauVSPARKGLRKLobattoIIIAB(s)` (2s); `TableauHPARKGLRK(1)`, `HSPARK(SPARKGLRK(s))`,
  `HSPARK(SPARKLob{ABC,ABD}(s))`, `TableauHSPARKLobattoIIIAIIIBpSymmetric(2)`.

## Findings

| # | Category | Method(s) | Observed | Root cause | Action |
|:--|:---------|:----------|:---------|:-----------|:-------|
| S1 | **B — bug (fixed)** | `TableauHSPARKGLRKLobattoIII{AB,BA,D,E}` | `BoundsError` (`s×s` matrix at `[i,σ+1]`) | `getTableauHSPARK` built the momentum projection coefficients `a_p_2`/`a_p_3` as `s×s` (`= g.a`) although the HSPARK-secondary residual indexes them over the `R = σ` projective stages; the position coefficients `a_q_2`/`a_q_3` were already the correct `s×σ`. **Incomplete port.** | **Fixed:** build `a_p_2`/`a_p_3` as the `s×σ` conjugate-symplectic partners of `α_q_2`/`α_q_3`, mirroring the `a_q_2`/`a_q_3` construction |
| S2 | **B — bug (fixed)** | all `HSPARKsecondary` (`TableauHSPARKLobattoIII*`, `TableauHSPARKGLRKLobattoIII*`) | `SingularException` | The null-vector residual/component code was commented out while `Cache` still allocates the `μ` unknown (`cache.jl:178`), leaving an unconstrained zero row/column in the Jacobian. The working `VSPARKsecondary` has the identical block active. | **Fixed:** re-enabled the null-vector `components!`/`residual!`/`initial_guess!` blocks (matching `VSPARKsecondary`) and corrected the null-vector guess guard from the nonexistent `:λ` field to `hasnullvector`. Moved the zero pivot off the null-vector row — see S3 |
| S3 | B — deeper (still broken) | all `HSPARKsecondary` | `SingularException` (now in the `ω` secondary-constraint block) | After S1/S2 the remaining singularity is a residual degeneracy in the `ω`-averaged secondary-constraint rows of this EXPERIMENTAL Hamiltonian method. The (incomplete) manuscripts do not give enough to reconstruct the intended coefficients with confidence, so no speculative numeric change was made. | Kept `@test_broken` with a precise, updated root cause |
| S4 | B — latent (fixed) | `HPARK`, `HSPARKprimary` | none observed (`P = 1`) | The δ-constraint residual zeroed row `R-1` (hard-coded) while accumulating into row `i` inside a `for i in R-P+1:R` loop; they coincide only when `P = 1` (all tested methods). | **Fixed:** zero the same row `i` that is accumulated (`integrators_hspark.jl`, `integrators_hspark_primary.jl`); behaviour-neutral for `P = 1`, correct for `P > 1` |
| S5 | **A — inherent** | `SPARKGLVPRK(2)`, `TableauHPARKGLRK(2)` | order 2 (docs 2s = 4) | `R(∞) = (-1)^{s+1} = -1` at `s = 2` violates the projection symplecticity conditions (prediction 2). This is the source's own `# maybe problem with R∞?` TODO. | `@test_broken` at order 4 |
| S6 | **A — inherent** | `SPARKLobattoIIIAIIIB(2/3/4)`, `SPARKGLRKLobattoIII{AIIIB,BIIIA}(s)` | order-reduced (`(3)`→~2.4, `(4)`→~3.5; GLRK-Lobatto→~1) | Symplectic Lobatto pair enforcing the constraint at the solution on the degenerate Lagrangian — reduced order (predictions 1 + the known symplectic-RK-on-degenerate-Lagrangian reduction). `SPARKLobattoIIIAIIIB(2)` fails the solve outright. | `@test_broken` / documented order-reduced tolerances |
| S7 | **A — inherent** | `SPARKLobattoIIIBIIIA(2/3/4)`, `TableauHPARKLobattoIII{AIIIB,BIIIA}(2/3)`, `TableauVSPARKLobattoIIIBIIIApSymmetric(3)` | divergence (`NonlinearSolverException`, or error 0.15–20 at `T = 0.1`) | Full divergence predicted by prediction 1; a different solver / line search / iteration cap does not help (confirmed here and in the third pass). | `@test_broken` |
| S8 | **A — inherent** | `VSPARK(SPARK{GLRK,LobABC,LobABD,LobattoIIIAIIIB,LobattoIIIBIIIA}(2))` | `SingularException` at `s = 2` | Degenerate stage system at the lowest stage count (prediction 3). `VSPARK(SPARKLobABD(2/3))` singular at `s = 2, 3`. | `@test_broken` |
| S9 | C — limitation (documented) | `SPARKMethod` (internal stages) | none on the degenerate test problem | `initial_guess!` zeroes the internal velocity `Vi` and `components!` never recomputes it (the `v̄` call is commented out), so `ϑ`/`f` are evaluated at `Vi = 0`. Harmless for degenerate Lagrangians (the `v`-term vanishes), but wrong in general — the code flags this in-comment. No non-degenerate DAE test problem exists here to validate a change, so a fix would be unverifiable. | Documented; not changed |

## Static consistency

The `components!`/`residual!`/`update!` stage equations match the docstrings and
the manuscript stage/update equations for every family; the `SLRK` and
`VSPARKsecondary`/`HSPARKsecondary` docstrings transcribe the paper's `(s,σ)`
constructions faithfully. The tableau `o` fields propagate correctly from the
underlying RungeKutta tableaus (`g.o`, `min(...)`), with the two closed-form orders
(`lobatto_gauss_coefficients` `o = s²` and `SLRKLobattoIII` `o = 2s-2`) matching
the code. The `docs/src/integrators/spark.md` page was expanded from a bare table
to a description of the method families, their DAE targets, the two Butcher-tableau
pairs, and the symplecticity/order caveats above.

## Fixes applied

* `src/spark/tableaus_hspark_secondary.jl` — build `a_p_2`/`a_p_3` as `s×σ`
  conjugate-symplectic matrices (finding S1); removed the stale "not used anymore"
  comment.
* `src/spark/integrators_hspark_secondary.jl` — re-enabled the null-vector
  `components!`/`residual!` blocks and fixed the `initial_guess!` guard (S2).
* `src/spark/integrators_hspark.jl`, `src/spark/integrators_hspark_primary.jl` —
  corrected the δ-constraint residual row index `R-1 → i` (S4).

## Tests

* **New** `test/verification/spark_convergence_tests.jl` (wired into
  `test/runtests.jl`): 64 order assertions + 11 `@test_broken` deficiencies.
* `test/spark/spark_integrators_tests.jl` — the HSPARK-secondary testset comment
  was updated to the S1–S3 root cause and the calls wrapped in `muffle` (the fixed
  methods now iterate in the solver before failing, rather than aborting
  immediately). Tally unchanged at **133 pass / 45 broken**.

The integrator machinery is correct: every SPARK family that is *not* subject to an
inherent instability (SLRK, SPARK-Gauß, SPARK-Lobatto ABC/ABD, VSPARK/VPARK
symplectic-projection, VSPARK-secondary, HSPARK-Gauß/Lobatto ABC/ABD) reaches
exactly its documented order. The non-converging cases are, with the exception of
the EXPERIMENTAL HSPARK-secondary family (S3), **inherent method properties**
(order reduction / divergence / singular stage systems predicted by the
manuscripts), not implementation defects.

---

# Fourth pass — DGVI, FLRK and PGLRK modernisation and verification

The three integrator families that the passes above recorded as modernisation blockers
have been ported to the method-based `GeometricIntegrator` architecture, registered as
methods, documented, and covered by the same static + dynamic verification standard as
the rest of the package. Nine genuine defects were found in the process, four of them
numerical.

## Summary of findings

| # | Severity | Method(s) | Issue | Root cause | Action |
|:--|:---------|:----------|:------|:-----------|:-------|
| 13 | correctness (**fixed**) | `CoefficientsPGLRK` | The constructor could not run at all — `TypeError` on every call | `new(...)` passed 11 arguments to a 14-field struct: `@CoefficientsRK` gained `â, b̂, ĉ` in `bb8c298b`, so `P` landed in `â` and the `Q` *matrix* in `b̂::Vector{T}` | **Fixed** by supplying zero `â/b̂/ĉ`, as `CoefficientsARK` already did |
| 14 | correctness (**fixed**) | `CoefficientsPGLRK` | `Legendre` was never imported into `Integrators` | `src/Integrators.jl` imported only `Basis` and `nbasis` from CompactBasisFunctions | **Fixed**; `nodes`/`weights` also qualified as `QuadratureRules.…` |
| 15 | doc bug (**fixed**) | `docs/src/audit.md` | Claimed `QuadratureRules` `nodes`/`weights` "now return matrices" | They return `Vector{Float64}`; the 0-based `leg_basis[nodes, j-1]` indexing was correct all along | **Corrected** above |
| 16 | correctness (**fixed**) | `FLRK` | The momentum update used the projection field `g = (∇ϑ)ᵀv` as the stage force instead of `f = (∇ϑ)ᵀv − ∇H`, so `ṗ` was missing `−∇H` | Legacy `integrate_diag_flrk!` called `equs[:g]` | **Fixed** to `equations(int).f`. Confirmed by reverting: `max|p − ϑ(q)|` degrades from 1E-15 to **9.9E-2** at every order, while `q` stays correct — which is why the old test, checking only `q`, never caught it |
| 17 | correctness (**fixed**) | `FLRK` | Would have been silently replaced by a plain `IRK` | `src/integrators/rk/abstract.jl:60` defines `initmethod(::RKMethod, …) = RK(method, TT)`; the historical `FLRK <: RKMethod` is caught by it | **Fixed** by making `FLRK <: LODEMethod`, which also drops the false `isodemethod`/`ispodemethod`/`ishodemethod` traits |
| 18 | correctness | `CoefficientsPGLRK(2)` | A two-stage method is inconsistent by construction | The skew perturbation occupies `W[2,1]`/`W[1,2]`, the one pair that coincides with the order-determining `ξ₁` entries of `X`. Verified: `bᵀA = 0.5` and `A·1 = 1.0` at `s = 2` versus ~1E-16 for `s ≥ 3`, and a slot-by-slot sweep shows `(2,1)` is the *only* offending slot at every `s` | Constructor now requires **`s ≥ 3`** with an explanatory message. This also explains the old disabled test's `6E-6` tolerance at `s = 2` against `2E-12` at `s = 3` |
| 19 | correctness | `VPRKpTableau` | Needs two distinct stage-count bounds | `D` multipliers occupy slots `(s,s−1) … (s−D+1,s−D)`: `s ≥ D+1` merely to fit into `W`, `s ≥ D+2` to avoid the `(2,1)` slot and keep full order | `s ≥ D+1` asserted; `s ≥ D+2` documented. Measured on the 2-dimensional problem: `s = 3` gives 2.1E-10 (worse than `VPRK(Gauss(3))`), `s = 4` gives 8.9E-16 |
| 20 | correctness (**fixed**) | `VPRKpTableau` | The legacy file referenced unbound names `method`, `solstep`, `problem` as if they were globals, and needed a non-dependency (`NLsolve`) | Bit-rotted mid-refactor | Rewritten as `src/integrators/vi/vprk_ptableau.jl` with the multipliers folded into a single coupled Newton system of size `D(s+1)`; `NLsolve` no longer needed |
| 21 | correctness | `DGVI`, `DGVIP0`, `DGVIP1` | Order capped at `2⌊s/2⌋` instead of `2s` | The trapezoidal flux is evaluated at the *nodal* value `qₙ`, which is only second-order accurate and limits the whole scheme. Consistent with the manuscript's own `% TODO The jump condition is not quite right` | Recorded as `@test_broken` against the nominal `2s`; the achieved orders are asserted |
| 22 | observation | the disabled DGVI tests | Could never have run | They used the *regular* harmonic oscillator, on which `DGVI`'s closure row collapses to `q − p` — independent of every unknown, hence a singular Jacobian. DGVIs need a **fully degenerate** Lagrangian | Tests moved to `LotkaVolterra2d.iodeproblem_dg` / `LotkaVolterra2dGauge` |
| 23 | dependency | `GeometricProblems` **0.6.24** (the pinned version) | `LotkaVolterra2d.iodeproblem_dg_gauge` fails `check_methods`; `PointVortices{,Linear}.lodeproblem_formal_lagrangian` fails to construct | The former's `g` closure has arity `(g,t,q,λ,params)` while `IODE` requires `(g,t,q,v,λ,params)`; the latter still passes `Ω=`/`∇H=` kwargs the current `LODEProblem` signature rejects | **Resolved.** Both were already fixed upstream in v0.7.0 (2026-07-12, `01e8026` "Standardize the problem interface across all example modules"), so this was never an upstream defect. The compat entry is now `GeometricProblems = "0.7"` — 0.6 is excluded outright, having carried serious bugs — which removes both workarounds: the DGVI gauge test uses `iodeproblem_dg_gauge(; κ)` directly, and the FLRK suite gains `PointVortices.lodeproblem_formal_lagrangian`, the problem written for that formulation |

| 24 | correctness (**fixed**) | `PGLRK` | `issymplectic` returned `true` | The *tableau* `a(λ)` is symplectic for every fixed `λ`, but `λ` is solved from the energy condition at every step, so the step map is that tableau composed with a state-dependent parameter choice and its Jacobian carries an extra rank-one term `(∂Ψ/∂λ)(∂λ/∂q)ᵀ`. Exact energy conservation and symplecticity are in any case mutually exclusive for a general Hamiltonian system unless the method reproduces the exact flow (Ge & Marsden, 1988), and `\|ΔH\|` is measured at 1.1E-15. Note the reference claims preservation of *quadratic invariants*, not symplecticity | **Fixed** to `missing`. Measured on a nonlinear pendulum at `Δt = 0.8` by finite-differencing the step map: the defect `\|JᵀΩJ − Ω\|_∞` converges to a constant **3.06E-7** as the FD step is refined (ε = 1E-3 … 1E-5), whereas `Gauss(3)`'s falls as `O(ε²)` to 7.3E-12 — i.e. Gauss's residual is FD noise and PGLRK's is real. At `Δt = 1.6` the two coincide, consistent with the λ solve finding no root there and falling back to plain Gauss |
| 25 | correctness (**fixed**) | `DGVIP0` | The documented `nbasis == nnodes` precondition was not enforced anywhere | The boundary one-forms are reconstructed by contracting `r∓` (length `S`) with the one-form sampled at the quadrature nodes (length `R`). The docstring stated the constructor required `S == R`, but no such check existed | **Fixed**: `check_basis_quadrature`, a per-variant hook on the generated DGVI constructors, asserts it. Previously `DGVIP0(Lagrange(nodes(GaussLegendreQuadrature(3))), GaussLegendreQuadrature(4))` constructed silently and then raised a bare `BoundsError` inside the first step for `R > S`, and would have truncated the reconstruction silently — a wrong answer with no error — for `R < S` |

Findings 13, 14, 16, 17, 20, 24 and 25 were implementation bugs in this repository and are
**fixed**. Findings 18, 19 and 21 are inherent properties of the methods, now asserted
or recorded. Finding 15 corrects this report; 22 corrects the test suite; 23 was a stale
dependency bound, now restricted to `GeometricProblems = "0.7"`.

Note that 24 is the one place where the port asserted a structural property that does not
hold: the same standard applied to `FLRK` (`issymplectic = missing`, because the reference
proves nothing) had been applied to `PGLRK`'s *tableau* rather than to the method.

## A structural result used throughout

The `CoefficientsPGLRK` family `a(λ) = a + λA` with `A = P W Q` has two independent
properties, both now asserted to machine precision in
`test/methods/pglrk_coefficients_tests.jl`:

* `a = P X Q` reproduces the Gauß tableau **exactly** (measured ≤ 2.4E-16 for `s = 3…6`),
  and `b`, `c` *are* the Gauß weights and nodes. This is a decisive check on the whole
  W-transformation construction.
* `B A` is **skew** for `B = diag(b)`, because `Pᵀ B P = I` (the Legendre basis being
  orthonormal at the Gauß nodes). Hence `a(λ)` satisfies the symplecticity condition
  `bᵢaᵢⱼ + bⱼaⱼᵢ = bᵢbⱼ` for **every** `λ` and every choice of nonzero slot in `W`.

This is *orthogonal* to the order question: the `(2,1)` slot breaks `bᵀA = 0` and
`A·1 = 0` — and hence consistency — while leaving symplecticity intact. A two-stage
PGLRK would therefore be a symplectic method of reduced order, which is why it is
rejected outright rather than silently allowed.

Both statements are about the **tableau at fixed `λ`**, and neither transfers to a method
that chooses `λ` from the state — see finding 24. `VPRKpTableau` is subject to the same
caveat, and already reports `issymplectic = missing`.

## Verified-correct behaviour (dynamic confirmation)

**FLRK** — the position phase is an implicit Runge-Kutta method applied to `q̇ = v̄(q)`
and reproduces `integrate(ode, Gauss(s))` to 1 ulp (bit-identically at fine `Δt`; the
1-ulp difference at coarse `Δt` comes from FLRK seeding Newton from `v̄` on top of the
Hermite extrapolation). Convergence order measured at **exactly `2s` in both `q` and
`p`** — 2.00/4.01/5.99 and 2.00/4.00/6.00 for `s = 1,2,3` — which is worth recording
because the reference proves no order theorem at all. Energy drifts slowly, as its own
modified equation predicts, so no energy assertion is made.

**PGLRK** — order matches the underlying Gauß method (5.69 / 7.74 against Gauß's
5.88 / 7.87 at `s = 3,4`; both sit slightly under nominal on this problem and step
range). Energy conservation is the point of the method and is confirmed on the
*nonlinear* Lotka–Volterra problem: **3.6E-15 against Gauß's 7.3E-10** at `s = 3`, five
orders of magnitude. The harmonic oscillator cannot distinguish them — Gauß already
conserves quadratic invariants exactly — so it is not used for that test. What is
conserved is the energy and *not* the symplectic form: the measured symplecticity defect
is 3.06E-7 at `Δt = 0.8`, of the same order as the energy error Gauß commits there
(4.1E-7), which is the trade the method makes (finding 24).

The `λ` solve costs roughly 8 additional nonlinear stage solves per step, so `PGLRK(3)` is
about 9× the cost of `Gauss(3)` (16.9 ms against 1.9 ms over 100 Lotka–Volterra steps).
Two of those were removed by not pre-probing the bracket endpoints that
`SimpleSolvers.bisection` evaluates itself. Replacing the bisection with a secant on `λ`
would cut most of the remainder and is recorded as possible future work; the residual is
smooth in `λ`, so it should converge in 3–5 iterations rather than one solve per bit.

**VPRKpTableau** — enforces the Dirac constraint `ϑ(qₙ) − pₙ = 0` to 4.4E-16 … 8.0E-15,
and recovers on a degenerate Lagrangian the order that a plain VPRK loses there: plain
`VPRK(Gauss(s))` shows the documented reduction to ≈ `s` (3.6 / 4.1 / 5.8 for
`s = 3,4,5`) while `VPRKpTableau` reaches 5.3 / 9.5 / 9.8. A clean order *number* is not
reproducible, however — the multiplier introduces a `Δt`-independent error floor
(≈1E-11 at `s = 4`, ≈1E-13 at `s = 5`) and the error sequence is non-monotone, so the
fitted slope swings with the step window. Tightening the solver tolerances does not move
the floor, so it is a property of the method. Recorded as `@test_broken`; the method has
no published order theorem.

**DGVI** — all five variants run and converge on `LotkaVolterra2d.iodeproblem_dg`. They
split sharply by flux discretisation (measured at `T = 1`, `Δt = 1/4 … 1/32`):

| `s` | `DGVI` | `DGVIP0` | `DGVIP1` | `DGVIPI` | `DGVIEXP` |
|:----|:-------|:---------|:---------|:---------|:----------|
| 2 | 1.94 | 2.01 | 1.94 | 1.97 | 2.01 |
| 3 | 2.01 | 2.01 | 2.01 | **5.95** | **5.97** |
| 4 | 3.74 | 3.86 | 3.74 | 5.97\* | 6.00\* |

(\* machine-precision limited.) The two variants that evaluate the one-form at an
*average* of the one-sided limits reach the full `2s`; the three that use the nodal value
are capped at `2⌊s/2⌋` (finding 21). At `s = 4` and `Δt = 0.05` that is 8.0E-14 and
5.2E-14 against 2.6E-8 — seven orders of magnitude. Note that `DGVIEXP`, the variant
carrying no derivation at all, is among the two best performers. `DGVI` and `DGVIP1`
agree to ~1E-12, as expected since `DGVIP1` only adds a continuity constraint and a
final projection on top of `DGVI`'s equations. The gauge-transformed problem
(`LotkaVolterra2dGauge`) is also verified.

## Fixes and changes applied

* `src/integrators/rk/pglrk_coefficients.jl` → **`coefficients_pglrk.jl`**, with the two
  constructor bugs fixed, `s ≥ 3` enforced, `nstages`/`eachstage`/`order`/`isapprox`
  added, a readable `show`, and an in-place `getTableauPGLRK` that accepts dual numbers.
  `CoefficientsPGLRK` and `getTableauPGLRK` are now exported.
* `src/integrators/rk/integrators_flrk.jl` — ported; the `ϑ/P/F/G/J/A` workspace moved
  from integrator fields into the cache (as integrator fields they defeated the
  dual-number `cache(int, ST)` mechanism entirely); Jacobians via
  `SimpleSolvers.Jacobian` rather than a direct ForwardDiff dependency (removed from
  `Project.toml` in `a46ef1a5` and no longer loadable); `LinearAlgebra.lu!`/`ldiv!` in
  place of the removed `LUSolver`.
* `src/integrators/rk/integrators_pglrk.jl` — ported. `λ`, the reference energy and the
  working tableau live in a mutable cache; the outer λ solve uses
  `SimpleSolvers.bisection` (which *does* exist as a scalar root finder, unexported)
  guarded against the no-sign-change case, where it would otherwise return an endpoint,
  i.e. the largest admissible perturbation. **Note for future work:** the working
  tableau must be built with `Matrix{ST}(coeff.a)`, not `convert(Matrix{ST}, coeff.a)` —
  the latter returns the argument itself when the type already matches, which aliased the
  cache to the method's own tableau so that every update accumulated into and permanently
  corrupted it. The symptom was PGLRK conserving energy *worse* than plain Gauß.
* `src/integrators/vi/vprk_ptableau.jl` — new, replacing the deleted legacy
  `src/integrators/vprk/integrators_vprk_ptableau.jl`.
* `src/integrators/dgvi/integrators_dgvi_common.jl` — new shared layer: one abstract
  `DGVIMethod`, one superset coefficient block (generated once for the five method
  types), one pruned superset cache, and the shared `initial_guess!`, stage
  computations, `update!` and `integrate_step!`. The five per-variant files now hold only
  a docstring, `description`, the flux, the residual and the state hooks. The cache
  dropped 13 of the legacy 41 `D`-vectors that were computed but never read, and the
  six-per-residual `zero(cache.q)` allocations inside the Newton loop became one
  permanent field.
* The carried-over jump values live in the **`DT` cache only** and are always reached as
  `cache(int).state`, never `cache(int, ST).state` — they are constant with respect to
  the solver unknowns and must not be dual-typed. `DGVI` needs no state at all: `p`
  carries the jump information, making it a genuine `(q,p)` map.
* Latent DGVI bugs removed in passing: `params.q⁻` was read but never written (harmless
  only because the fields depending on it were unread); `DGVIP0` used `q̄⁺`/`Θ̅⁺` that were
  never computed; `DGVIP0` looped `for i in 1:S` over an array of length `R`, a
  `BoundsError` whenever `S > R`; and the experimental variants accumulated `t` by
  `params.t += Δt` from zero rather than taking it from the solution.
* `docs/legacy.md` (untracked) deleted, its one unique row — `VPRKpTableau` — folded into
  `docs/src/integrators/vprk.md`. The DGVI documentation page is re-enabled in
  `docs/make.jl`, two notational bugs in it corrected (the `r∓` reconstruction indices,
  which contradicted the same page's own variation formulae, and the missing transpose in
  `∇ϑ`), and a section added mapping the five variants onto the framework.
  `docs/src/modules/integrators.md` gained the new files and lost four `Pages` entries
  naming files that do not exist.
* `test/methods/test_list.jl` (untracked, never run, referenced from nowhere) replaced by
  `test/methods/method_list_tests.jl`, a real testset that sweeps every registered method
  through `order`/`isexplicit`/`is*method`. Registering a method without its trait
  overloads breaks the *docs build* rather than the test suite, so this guards a gap that
  nothing else covered. `test/integrators/test_show.jl`, likewise orphaned, was repaired
  (it still passed a raw `Tableau` to `GeometricIntegrator`) and wired into `runtests.jl`.

## Tests

New: `test/methods/pglrk_coefficients_tests.jl` (102 assertions),
`test/methods/method_list_tests.jl` (135), `test/verification/flrk_convergence_tests.jl`
(21), `test/verification/pglrk_convergence_tests.jl` (15 + 1 broken),
`test/verification/dgvi_convergence_tests.jl` (20 + 3 broken).

Re-enabled with empirically measured tolerances: the FLRK block in
`rk_implicit_integrators_tests.jl`, the PGLRK block in `rk_integrators_tests.jl`, all
five DGVI variants in `galerkin_integrators_tests.jl` (35 assertions, on a degenerate
problem), the `VPRKpTableau` block in `projections_vprk_tests.jl`, the
`CoefficientsPGLRK` check in `spark_tableaus_tests.jl`, and the FLRK entries in
`methods/runge_kutta_methods_tests.jl` and `integrators/test_show.jl`.

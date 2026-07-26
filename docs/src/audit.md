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
the manuscript stage/update equations for every family; the
`VSPARKsecondary`/`HSPARKsecondary` docstrings transcribe the paper's `(s,σ)`
constructions faithfully. **Correction (fifth pass):** the claim originally made here
that the `SLRK` docstring does so too was wrong — it was a verbatim copy of the
`VSPARKsecondary` docstring and described a different scheme. See finding S11 below.
The tableau `o` fields propagate correctly from the
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

---

# Fifth pass — `SLRK` against the SLRK manuscript

The fourth pass verified the SPARK submodule against the two *SPARK Methods for …*
drafts. `SLRK`, however, has its own source manuscript —
`/Users/mkraus/Datashare/Papers/Symplectic Lobatto Runge-Kutta Methods for Degenerate
Lagrangian Systems/slrk_degenerate_lagrangians.tex` (Dec 2019) — which the fourth pass
did not use. This pass checks `src/spark/integrators_slrk.jl` and
`src/spark/tableaus_slrk.jl` line by line against it.

The manuscript is an incomplete working draft: Sections 1–3 are largely empty, the
`.bib` is 0 bytes, six `\eqref`s point at labels that are never defined, and the
author left a `% TODO Check signs!`. It was judged critically rather than taken as
ground truth — and the most consequential finding is a gap in *its* proof, not in the
code.

## What matches

Notation map (code → manuscript): `q.a → a¹`, `q̃.a → a²`, `p.a → ā¹`, `p̃.a → ā²`,
`q.b → b¹`, `q̃.b → b²`, `p.b → b̄¹`, `p̃.b → b̄²`. All `SLRKLobattoIII*` set `q̃ = q`,
`p̃ = p`.

| Manuscript | Code | Verdict |
|:---|:---|:---|
| Stage equations for `P`, `Q` | `residual!` rows `+1`/`+2` | match |
| Update rules for `q`, `p` | `components!`, `update!` | match |
| `F = −∇H + ∇ϑ·V`, `G = ∇ϑ·Λ` | `equations.f`, `equations.g` | match |
| Primary constraint at **every** stage, un-averaged (the paper's central idea) | row `+3` | match |
| `0 = Σⱼ ωᵢⱼ ψⱼ + ω_{i,s+1} ϕ(qₙ₊₁,pₙ₊₁)`, `i = 1..s`, `ω ∈ ℝ^{s×(s+1)}` | rows `+4` | match |
| `0 = Σᵢ dᵢ Vᵢ` | null-vector block | match |
| `ψ = F̃ − Ṽ·∇ϑ(Q̃)` | `equations.ψ(…, Vp, Fp)` | match |

Verified numerically (`scripts/slrk_verification.jl`, steps 1a/1b):

* **Symplecticity conditions of the tableaus.** Every `X̄` in RungeKutta.jl is
  `symplectic_conjugate_coefficients(a,b) = bⱼ/bᵢ (bᵢ − aⱼᵢ)`, so `b̄ᵢaᵢⱼ + bⱼāⱼᵢ = b̄ᵢbⱼ`
  and `b = b̄` hold by construction. Worst deviation over all six constructors and
  `s = 2,3,4`: **2.8e-17**. With `q̃ = q`, `p̃ = p` the manuscript's eight conditions
  collapse to these two.
* **The `ω` matrix.** `lobatto_ω_matrix(s)` builds `ω = L⁻¹R`; since `L` is invertible,
  `ωX = 0 ⟺ RX = 0`, i.e. exactly the `s−1` Lobatto-IIIA-averaged secondary constraints
  plus `ϕ(qₙ₊₁,pₙ₊₁) = 0`. Confirmed by cross-nullspace residual ≤ 1e-15 for
  `s = 2,3,4`; for `s = 2`, `ω = [1 1 −1; 0 0 1]`.
* **Primary-constraint preservation**, the property the family exists for:
  `max|ϕ(qₙ,pₙ)| ≈ 1.8e-15` over 1000 steps at `Δt = 0.1`, for every constructor and
  every `s`.

## Findings

| # | Category | Issue | Root cause | Action |
|:--|:---------|:------|:-----------|:-------|
| S10 | **B — bug (fixed)** | The null-vector multiplier `μ` was added to **two** residual rows: the momentum-stage row (`+2`, Z-space) *and* the primary-constraint row (`+3`, P-space). Since `P = p + h·Z`, the same coefficient carries an extra factor `h` in the first, so the two combine to `p + h(…) − ϑ(Qᵢ) = (1−h)·μ·dᵢ/b̄ᵢ`. The stage Jacobian is therefore **exactly singular at `Δt = 1`** and ill-conditioned near it. Present since the original commit `cb1d6406` (2020-02-21). | A copy-paste merge of the two conventions used elsewhere in the submodule: `vpark`/`vspark`/`spark` put `μ` on the `Φ` row only, `vspark_secondary`/`hspark_secondary` on the `Z` row only. `SLRK` was the only integrator doing both. | **Fixed:** dropped the `Z`-row term, kept the `Φ`-row term |
| S11 | **B — doc bug (fixed)** | The `SLRK` docstring was a **verbatim copy** of the `VSPARKsecondary` docstring and described a different method: two stage sets (`s` internal + `σ` projective), three vector fields `F¹ = ∇H`, `F² = ∇ϑ·V`, `F̃³ = ∇ϕ·Λ`, three momentum weight sets. `SLRK` has one stage set of size `s`, `F = ∂L/∂q`, `G = ∇ϑᵀΛ` and two weight sets, and the `μ`/`d` equations were undocumented. Inherited errors: `F¹ = ∇H` missing its minus sign, `Σⱼ ωᵢⱼ Ψ̃_{n,i}` should be `Ψ̃_{n,j}`, range `i = 1..σ−1` instead of `1..s`, `ω_{i,s+1} ϕ(qₙ₊₁,pₙ₊₁)` term omitted. | Copy-paste; the fourth pass asserted the opposite (corrected above). | **Fixed:** docstring rewritten from the code and the manuscript |
| S12 | **C — API trap (documented)** | `SLRK` requires the `LDAEProblem`'s `f` to be the **full** `∂L/∂q`, because `F` feeds directly into both the `p` update and `ψ`. `VSPARKsecondary` requires the *opposite* (`f = −∇H`, reconstructing `∂L/∂q = f + g(V)`). Passing a `VSPARKsecondary`-style problem to `SLRK` silently drops the `∇ϑ·v` term. Nothing said so; the docstring said the reverse. | Two conventions coexisting in one submodule — hence the two constructors `ldaeproblem_slrk` / `ldaeproblem` in `GeometricProblems.LotkaVolterra2d`. | **Documented** prominently in the docstring and `docs/src/integrators/spark.md` |
| S13 | **B — bug (fixed)** | `SLRKLobattoIIIC̄C` was registered under `Symbol("SLRKLobattoIIICIIIC̄")` — the same name as `SLRKLobattoIIICC̄`. | Copy-paste. | **Fixed** to `Symbol("SLRKLobattoIIIC̄IIIC")` |
| S14 | C — robustness (fixed) | `components!` evaluated `ϕ` with `C.ṽ`, which is never assigned (relies on the cache staying zero); `components!` used `p.c[i]` for the stage time while `initial_guess!` used `p̃.c[i]`, although `Qᵢ` is generated by the `q` tableau. Behaviour-neutral today (all Lobatto pairs share nodes, and `ϕ = p − ϑ(q)` ignores `v`). | Drift. | **Fixed:** both use `q.c`, `C.ṽ` zeroed explicitly, both commented |
| S15 | **A — inherent (manuscript proof is wrong)** | **`SLRK` is not symplectic.** The Poincaré invariant `∮p·dq` drifts secularly. | See below. | Documented in the docstring, in `docs/src/integrators/spark.md`, and corrected in the manuscript |

### S15 in detail

The manuscript's Theorem claims exact preservation of the noncanonical symplectic
two-form. Its proof does not close. Reconstructing it in the code's notation, the
first-order terms reduce to

```
T₁+T₃ = h Σⱼ (bⱼ − b̄ⱼ) dϑ∧dVⱼ + h Σⱼ (bⱼ dⱼ/b̄ⱼ) dμ∧dVⱼ
T₂+T₄ = h Σⱼ (βⱼ − β̄ⱼ) dϑ∧dΛⱼ + h Σⱼ (βⱼ dⱼ/b̄ⱼ) dμ∧dΛⱼ
```

The first summand of each vanishes by `b = b̄`. The `dμ∧dV` term vanishes by the
null-vector constraint `Σᵢ dᵢVᵢ = 0`. **The `dμ∧dΛ` term has no counterpart** — nothing
constrains `Σᵢ dᵢΛᵢ`. In the manuscript this is exactly the point where the proof
defers to `\eqref{eq:projected_lobatto_symplecticity_conditions_constraint}`, *a label
that is never defined anywhere in the file*.

Measured (`scripts/slrk_verification.jl`, steps 4–7):

* **Control.** The same `∮p·dq` harness (spectral differentiation of the loop, 200
  points, radius 0.1) applied to genuinely symplectic variational integrators —
  `VPRKGauss(2)`, `VPRKLobattoIIIAIIIĀ(2)`, `VPRKLobattoIIIAIIIĀ(3)` on the `lode` form
  of the same problem — gives `|ΔI|/|I| ≈ 1e-14` after 1, 10 **and** 100 steps. The
  harness is clean.
* **`SLRK`.** `|ΔI|/|I|` is nonzero and accumulates linearly with the step count. For
  `SLRKLobattoIIIAB(2)` at `Δt = 0.1`: `2.6e-4` (1 step) → `2.0e-3` (10) → `1.6e-1`
  (100). Per-step it scales as `O(h^{p+1})`:

  | method | h=0.2 | h=0.1 | h=0.05 | h=0.025 |
  |:---|---:|---:|---:|---:|
  | `SLRKLobattoIIIAB(2)` | 1.5e-3 | 2.6e-4 | 3.6e-5 | 4.7e-6 |
  | `SLRKLobattoIIIAB(3)` | 2.2e-6 | 7.8e-8 | 2.6e-9 | 8.3e-11 |
  | `SLRKLobattoIIID(3)`  | 9.4e-8 | 8.6e-10 | 7.3e-12 | 1.9e-13 |

  An independent finite-difference measurement of `‖JᵀΩ(qₙ₊₁)J − Ω(qₙ)‖` on the
  constraint manifold agrees to within a few percent and is flat in the FD step `ε`
  over four decades, ruling out quadrature and solver noise.
* **Diagnosis confirmed.** Reading the stage variables out of the cache:
  `Σᵢ dᵢVᵢ ≈ 1e-16` (imposed) but `Σᵢ dᵢΛᵢ = O(h)` (uncontrolled). The product
  `h·|μ|·|Σd·Λ|` — the size of the surviving term — predicts the measured defect within
  ~10 % for most pairs (e.g. `SLRKLobattoIIIAB(2)` at `h = 0.1`: 2.7e-4 predicted vs
  2.57e-4 measured; `SLRKLobattoIIIBA(2)`: 6.0e-6 vs 6.5e-6).
* **Not caused by S10.** After the S10 fix the drift is unchanged digit for digit.

`SLRKLobattoIIID` and `SLRKLobattoIIIE` have the smallest defect. The methods remain
useful: they attain their documented order `2s−2` and preserve the primary constraint
exactly, which is what variational integrators cannot do. They are *approximately*
symplectic, with the violation at the level of the local truncation error.

Note that this is consistent with the fourth pass's prediction 1 (symplectic *and*
constraint-at-the-solution is impossible): `SLRK` does enforce `ϕ(qₙ₊₁,pₙ₊₁) = 0`
through the last column of `ω`.

## Verification of the S10 fix

| check | before | after |
|:---|:---|:---|
| `cond(J)` at `Δt = 0.9 / 0.99 / 0.999 / 1.0`, `SLRKLobattoIIIAB(2)` | 1.5e2 / 2.2e3 / 2.3e4 / **1.6e11** | 3.0e1 / 4.5e1 / 4.6e1 / **4.7e1** |
| `integrate` at `Δt = 1` | `SingularException` | `NonlinearSolverException` (Newton divergence at a step size comparable to the period — no longer a linear-algebra failure) |
| `q(T)`, `p(T)` at `Δt = 0.1`, all six constructors, `s = 2,3` | — | unchanged to ≤ 2 ulp |

The fix is behaviour-neutral at usable step sizes, as predicted: for problems whose
`g`, `ψ` and `ϕ` do not depend on `p` beyond `ϕ = p − ϑ(q)`, the two variants differ
only by the reparametrisation `μ ↦ (1−h)μ`. It *does* change the reported internal
stages `Pᵢ` and would change results for any `ψ` or `g` that depends on `p`.

## Fixes applied

* `src/spark/integrators_slrk.jl` — dropped the duplicated `Z`-row null-vector term
  (S10) with a comment recording the convention; rewrote the docstring (S11, S12);
  aligned the stage-time tableau and zeroed `C.ṽ` explicitly (S14).
* `src/spark/tableaus_slrk.jl` — corrected the duplicated tableau name (S13).
* `docs/src/integrators/spark.md` — added `IntegratorSLRK` to the family table and a
  full "Specialised Lobatto Runge–Kutta methods (`SLRK`)" section: what distinguishes
  it from Jay's SPARK methods, the complete scheme, the `ω` construction, the null
  vector, the table of the six constructors, the `f = ∂L/∂q` warning with a usage
  example, and what the family does and does not preserve.
* `scripts/slrk_verification.jl` — the seven-step verification script used above.

## Manuscript corrections

Applied to `slrk_degenerate_lagrangians.tex` (the paper is under git; the file still
compiles, 10 pages, with no new errors or undefined references beyond the pre-existing
ones):

| # | Correction |
|:--|:-----------|
| M1 | The Remark's `Ω(Q̃)·Ṽ = ∇H(Q̃)` was wrong by a sign given `Ω := ∇ϑ − ∇ϑᵀ` — this is the author's own `% TODO Check signs!`. `Ω` redefined as `∇ϑᵀ − ∇ϑ`, with the derivation from `ψ = 0` written out and the TODO removed. Code and `GeometricProblems`' `ω` already used the correct convention. |
| M2 | Definition and proof contradicted each other over where `μ` lives — the definition put it in the `P` stage equation and demanded `ϕ(Qᵢ,Pᵢ) = 0` exactly, while the proof substitutes `dP̃ᵢ = dϑ(Q̃ᵢ) + (dᵢ/b⁴ᵢ)dμ`. Adopted the proof's version: `μ` moved into the constraint. |
| M3 | Proof gaps. (i) After the `Ṽ` block, "vanishes if `b² = b⁴`" only kills two of three summands; the third needs constraint (3c) — now derived, and its sign corrected. (ii) The `Λ̄` block deferred to a **label that is never defined**; the surviving term is now written out explicitly, the missing condition stated as (3c)'s `Λ` analogue, and the numerical evidence for the drift summarised. |
| M4 | The proof uses three stage sets, five weight vectors `b¹…b⁵` and three coefficient families that the Definition never introduces. Added an explicit notation map and showed the eight conditions collapse to the two standard partitioned ones. |
| M5 | The `μ` term was the only stage term without an `h`; noted that `μ = O(h)` and that scaling by `h` merely reparametrises it. |
| M6 | Constraint (3c) was written with `Ṽ_{n,i}`, undefined in the Definition → `V_{n,i}`. Also `ω_{is+1}` → `ω_{i,s+1}`. |
| M7 | `ω` denoted two different matrices. The `(s−1)×s` Lobatto-IIIA block renamed `ω^A`, and the `s×(s+1)` construction the code uses added as a Remark, with the `s = 2` example. |
| M8 | The empty `\cite{}`, the empty `.bib` and the six dangling `\eqref`s are now listed as explicit submission blockers in the (empty) Introduction rather than failing silently. |
| M9 | §5.2's Lagrangian was missing the `+ q₁q̇₂` term, which is why its equations of motion carried a spurious `1/(1+q₁q₂)` factor and did not match `GeometricProblems.LotkaVolterra2d`. Term restored, equations of motion corrected, and a Remark added explaining that `(0,q₁)` is not a gradient (so the truncated version is a genuinely different symplectic structure, not a gauge choice) and that the paper's `b`-sign convention is the negative of the package's. |
| — | The abstract's claim of exact symplecticity was removed and replaced with the constraint-preservation claim plus a flagged caveat. The Theorem now states the missing condition explicitly. |

## Tests

No regression. `spark_tableaus_tests.jl` 130 pass; `spark_integrators_tests.jl`
**133 pass / 45 broken**; `spark_convergence_tests.jl` **64 pass / 11 broken** —
identical to the fourth-pass tally.

## Fifth pass, addendum — the same verification on `LotkaVolterra2dSingular`

`GeometricProblems.LotkaVolterra2dSingular` has `ϑ = (log q₂/q₁, 0)`, which differs
from `LotkaVolterra2d`'s `ϑ = (q₂ + log q₂/q₁, q₁)` by the **exact** form `d(q₁q₂)`.
The two therefore describe the *same* continuous system — same `Ω = ∇ϑᵀ − ∇ϑ`, same
Euler–Lagrange equations, same `∮p·dq` on any closed loop — but give different
discrete methods. Repeating steps 2, 4, 6 and 7 on both separates what is a property
of the *system* from what is a property of the *discretisation*.

Confirmed first: the loop integral is gauge invariant, `∮p·dq = −0.031573928398` for
both, as it must be since `∮∇F·dq = 0`.

### The methods split into two groups

The split is predicted exactly by the ranks of the two coefficient matrices:

| method | `rank(a¹)` = `rank(q.a)` | `rank(ā¹)` = `rank(p.a)` | both singular | on the singular gauge |
|:---|:---|:---|:---:|:---|
| `SLRKLobattoIIIAB(s)` | `s−1` | `s−1` | **yes** | `SingularException` |
| `SLRKLobattoIIIBA(s)` | `s−1` | `s−1` | **yes** | singular / divergent |
| `SLRKLobattoIIICC̄(s)` | `s` | `s−1` | no | fine |
| `SLRKLobattoIIIC̄C(s)` | `s−1` | `s` | no | fine |
| `SLRKLobattoIIID(s)` | `s` | `s` | no | fine |
| `SLRKLobattoIIIE(s)` | `s` | `s` | no | fine |

Checked for `s = 2,3,4`. The two Lobatto IIIA–IIIB pairs are exactly the ones where
`a¹` *and* `ā¹` are both rank deficient — `IIIA` has a zero first row, `IIIB` a zero
last column, and each conjugate partner inherits the defect. On the standard gauge
this is harmless because `∇ϑᵀ` has full rank 2 and the missing direction is supplied
by the projection force `G = ∇ϑᵀΛ`. On the singular gauge `ϑ₂ ≡ 0`, so the second
column of `∇ϑᵀ` vanishes, `∇ϑᵀ` drops to rank 1, `G` can no longer reach the second
momentum component, and the stage system loses rank outright.

The failure is *not* the `Δt = 1` conditioning problem of S10 — it is present at every
step size:

| `cond(J)`, singular gauge | Δt=0.1 | Δt=0.5 | Δt=0.9 | Δt=1.0 |
|:---|---:|---:|---:|---:|
| `SLRKLobattoIIIAB(2)` | 3.1e12 | 1.0e17 | 2.7e17 | 9.3e11 |
| `SLRKLobattoIIIBA(3)` | 1.8e14 | 1.1e11 | 1.8e10 | 3.8e10 |

against 8.2e1 … 4.7e1 for the same methods on the standard gauge. `SLRKLobattoIIIBA(3)`
does complete a single step (agreeing with the standard gauge to 7 digits) but diverges
over ten: `‖JᵀΩJ − Ω‖` jumps to 1.1e2 … 3.0e0 with no convergence order, and the
`T = 100` run produces `NaN`.

This is the fourth pass's prediction 1 (symplectic Lobatto `IIIA-IIIB` pairs on a
degenerate Lagrangian) reappearing — now shown to be **gauge dependent**: the same
method is usable in one gauge and unusable in a gauge-equivalent one.

### For the four methods that work on both: no difference whatsoever

`CC̄`, `C̄C`, `D` and `E` are **gauge invariant to round-off**. Over `T = 1` at
`Δt = 0.1`, `max|q_standard − q_singular| ≤ 6.7e-16` for every one of them at `s = 2,3`.
Every conservation diagnostic follows suit — the numbers below are identical in both
gauges, digit for digit:

| diagnostic (Δt = 0.1) | `CC̄(2)` | `CC̄(3)` | `C̄C(2)` | `D(2)` | `D(3)` | `E(2)` |
|:---|---:|---:|---:|---:|---:|---:|
| `∮p·dq` drift, 1 step | 8.52e-05 | 1.72e-08 | 1.10e-04 | 4.38e-08 | 8.57e-10 | 1.50e-08 |
| `∮p·dq` drift, 100 steps | 6.90e-02 | 2.74e-05 | 5.70e-02 | 1.72e-03 | 5.74e-07 | 2.40e-04 |
| `‖JᵀΩJ − Ω‖`, 1 step | 8.63e-05 | 1.82e-08 | 1.11e-04 | 1.09e-07 | 8.41e-10 | 1.37e-08 |
| `max|ΔH|`, `T = 100` | 1.699e-01 | 6.816e-05 | 7.798e-02 | 9.052e-03 | 7.238e-06 | 7.909e-04 |
| `max|ϕ|`, `T = 100` | 1.3e-15 | 1.3e-15 | 1.6e-15 | 1.6e-15 | 1.3e-15 | 1.3e-15 |

This is not a coincidence: for these methods the discrete `q`-trajectory is identical,
and `∮p·dq = ∮ϑ(q)·dq` changes by `∮∇F·dq = 0` under a gauge transformation, so both
the invariant and its drift must agree exactly.

**Conclusion on the question of gauge dependence.** The S15 symplecticity defect is a
property of the *discretisation scheme*, not of the representation of the Lagrangian:
it is unchanged by a gauge transformation. What the gauge *does* change is
solvability — the two Lobatto IIIA–IIIB pairs break down on the singular gauge, so the
choice of one-form is not free in practice. `SLRKLobattoIIID` and `SLRKLobattoIIIE`
are the safest choices on both counts: full-rank coefficient matrices in both blocks,
and the smallest symplecticity defect.

An incidental observation: `Σᵢ dᵢΛᵢ` — the uncontrolled quantity behind S15 — is itself
gauge invariant for these four methods (e.g. `CC̄(2)` at `h = 0.2`: 1.70e-01 in both
gauges), while `|μ|` is not (3.5e-02 vs 1.8e-02). The product `h|μ||Σd·Λ|` is therefore
only an order-of-magnitude proxy for the surviving wedge term, as noted above; the
wedge itself is gauge invariant.

## Fifth pass, addendum 2 — can SLRK be made exactly symplectic?

The S15 defect is structural. This section records the attempt to remove it, a
no-go argument for the SLRK ansatz, and the verification of the alternatives.

### The obstruction

Write `a = a¹ = q.a`, `α = a² = q̃.a`, `ā = ā¹ = p.a`, `ᾱ = ā² = p̃.a` and
`b, β, b̄, β̄` for the corresponding weights. The `h²` terms of the wedge product give

```
(C1) b̄ᵢbⱼ = b̄ᵢaᵢⱼ + bⱼāⱼᵢ      (C2) b̄ᵢβⱼ = b̄ᵢαᵢⱼ + βⱼāⱼᵢ
(C3) β̄ᵢbⱼ = β̄ᵢaᵢⱼ + bⱼᾱⱼᵢ      (C4) β̄ᵢβⱼ = β̄ᵢαᵢⱼ + βⱼᾱⱼᵢ
```

and the `h¹` terms give `(C5) b = b̄`, `(C6) β = β̄`, plus the two multiplier terms.

1. **`(C1)`–`(C3)` force `αᵢⱼ = (βⱼ/bⱼ)aᵢⱼ` and `ᾱᵢⱼ = (βⱼ/bⱼ)āᵢⱼ`.** Hence `Q`, `P`,
   `q_{n+1}` and `p_{n+1}` depend on `V` and `Λ` **only** through
   `Wⱼ = Vⱼ + (βⱼ/bⱼ)Λⱼ`; only `Ψ` sees `V` separately.
   *Verified:* replacing `Σdᵢ Vᵢ = 0` by `Σdᵢ Λᵢ = 0` gives bit-identical trajectories,
   drifts and orders for all six families — two points on the same gauge orbit.
2. **The multiplier is forced.** `ϕ = 0` at `s` stages plus `ϕ(q_{n+1},p_{n+1}) = 0` is
   `(s+1)D` conditions on the `sD` unknowns `W` — over-determined by exactly `D`.
3. **The two leftover terms combine into `h·dμ ∧ Σⱼ dⱼ dWⱼ`.** Since `W` is fully
   determined by the remaining equations, `Σd·W` is a nontrivial function of the
   initial data, so the term vanishes only if `dμ ≡ 0`.
4. **`μ ≡ 0` is impossible.** It would require the endpoint constraint to be implied by
   the stage constraints, i.e. `a` and `ā = conj(a)` both stiffly accurate. But
   `ā` stiffly accurate ⟺ `a` has a zero last column ⟹ `a_ss = 0`, while `a` stiffly
   accurate ⟹ `a_ss = b_s`. So `b_s = 0` — impossible for Lobatto.
5. **Decoupling `Λ` is impossible.** `(C6)` is `β = β̄`, so letting `Λ` project only the
   position (`β̄ = 0`) forces `β = 0` and then `α = 0`: `Λ` disappears from the scheme.

Four modifications were implemented and measured (via `residual!` overrides; the repo
source was not changed):

| modification | result |
|:---|:---|
| drop `μ` and `d` entirely | **singular** for all six families, including `IIID`/`IIIE` where both tableau blocks are full rank |
| replace `Σdᵢ Vᵢ = 0` by `Σdᵢ Wᵢ = 0` | **singular** — `W` is already determined, the condition is redundant |
| replace `Σdᵢ Vᵢ = 0` by `Σdᵢ Λᵢ = 0` | identical results to the current scheme (point 1 above) |
| drop the endpoint constraint (`ω`'s last column) **and** `μ` | see below |

The last one is the informative case. Dropping the endpoint constraint removes the
over-determination, `μ` disappears, and the scheme becomes **exactly symplectic**:

| Δt = 0.1 | canonical ∮ 1 / 100 | noncanonical ∮ 1 / 100 | max\|ϕ\| 1 / T=100 | order |
|:---|:---|:---|:---|---:|
| `IIID(2)` current | 4.4e-08 / 1.7e-03 | identical | 1.1e-15 / 1.8e-15 | 2 |
| `IIID(2)` no endpoint | **2.2e-16 / 2.4e-15** | 5.5e-04 / 3.3e-02 | 3.3e-04 / 1.1e-01 | 2 |
| `IIID(3)` no endpoint | **2.2e-16 / 2.2e-15** | 1.2e-05 / 1.8e-07 | 5.7e-06 / 1.8e-05 | 4 |
| `IIIAB(2)` no endpoint | 6.6e-16 / 6.3e-14 | 0.0 / 9.0e-15 | 1.2e-01 / 1.0e+02 | **0** |

— but the primary constraint is then lost, the map leaves `{ϕ=0}`, and the
*noncanonical* invariant drifts instead. `IIIAB(3)`, `IIIBA`, `IIICC̄`, `IIIC̄C` go
singular without the endpoint constraint; `IIIAB(2)` loses consistency entirely.

**Conclusion (S16).** Within the SLRK ansatz the three properties — exact symplecticity,
primary constraint at every internal stage, and the constraint at the solution — cannot
be had together. This is the SLRK manuscript's own Section-3 remark (which it invokes to
criticise Jay's SPARK methods) applying to its own scheme.

### Long-time behaviour: secular vs bounded

`∮ϑ(q)·dq`, relative drift, `h = 0.1`, out to 3000 steps:

| | 1 | 10 | 100 | 1000 | 3000 | growth |
|:---|---:|---:|---:|---:|---:|:---|
| `VPRKGauss(3)` | 1.4e-06 | 1.4e-06 | 5.8e-08 | 7.0e-06 | 7.8e-06 | **bounded** (dips ×24 at 100) |
| `VPRKGauss(2)` | 1.8e-04 | 1.9e-03 | 4.2e-03 | 3.5e-03 | 3.5e-02 | bounded / non-monotone |
| `SLRKLobattoIIID(3)` | 8.6e-10 | 9.5e-10 | 5.7e-07 | 6.1e-06 | 1.9e-05 | **linear** (×10.6, ×3.04) |
| `SLRKLobattoIIIAB(3)` | 7.8e-08 | 2.5e-07 | 4.1e-04 | 4.3e-03 | 1.3e-02 | **linear** (×10.5, ×3.03) |
| `SLRKLobattoIIIAB(2)` | 2.6e-04 | 2.0e-03 | 1.6e-01 | 7.2e-01 | 9.2e-01 | monotone, saturating at O(1) |

VPRK's canonical invariant stays at round-off throughout (2.2e-16 → 6.3e-14 over 3000
steps) while its noncanonical one merely oscillates with `max|ϕ|`; SLRK's grows linearly
in the step count to three significant figures over two decades. Different mechanisms:
a symplectic map measured with a non-invariant quantity, versus genuine non-conservation.

### S17 — the VSPARK projection methods are neither provably nor empirically symplectic

The natural alternative to SLRK is the projective-SPARK family (a symplectic inner method
plus a projection). Its theorem (*SPARK Methods for Degenerate Lagrangian Systems*,
`sec:projective_spark`) requires four tableau conditions, `b³ = b¹`, `b⁴ = b²`, **and**
`P̃_{n,i} = ϑ(Q̃_{n,i})` at *every* projective stage. That manuscript states plainly:

> "We find that none of the constructions satisfies all symplecticity conditions."

Confirmed, and it is not always `cond4` and not always `R(∞)`. Over the ten Gauß-inner
projections at `s = 1,2,3` (`scripts/vspark_projection_symplecticity.jl`, step 1):

| construction | condition that fails |
|:---|:---|
| `pSymplectic` | `cond4` at odd `s` (`R(∞) = -1`); holds at `s = 2` |
| `pSymmetric` | `cond4` = 0.25 at **every** `s` — its own docstring says so |
| `p{Modified,}Midpoint` | `cond2`/`cond3` at `s = 1,3`; hold at `s = 2` |
| `pInternal` | `cond2`/`cond3` at **every** `s` |
| `pModifiedInternal` | `cond2`/`cond3` at `s = 1,3` |
| `p{Modified,}LobattoIII{AIIIB,BIIIA}` | **all four hold at every `s`** |

For the last group the *only* remaining gap is the per-projective-stage requirement.

#### The earlier "empirically exactly symplectic" claim was a test-problem artefact

An earlier revision of this section recorded that the Gauß-based projections behave as if
exactly symplectic, on the strength of a one-step defect flat at round-off across a factor
of ten in `h`. That was measured on `LotkaVolterra2d` alone. Extending to a second
degenerate system — the massless charged particle, in both gauges — **overturns it.**

One-step canonical `∮p·dq` defect, `h = 0.5 → 0.05`, `N = 128` loop points:

| method | LotkaVolterra2d | LV2dSingular | MasslessChargedParticle | MCPSingular |
|:---|:---|:---|:---|:---|
| `GLRK(1)pSymplectic` | round-off | **fails** | **1.6e-05 → 1.0e-09**, `O(h⁴)` | **fails** |
| `GLRK(2)pSymplectic` | round-off | 5.6e-06 → r.o. | **3.7e-10 → 4.6e-15** | round-off |
| `GLRK(3)pSymplectic` | round-off | **2.7e-03**, `O(h⁴)` | **6.5e-11 → 1.8e-15** | **1.4e-05**, `O(h⁷)` |
| `GLRK(1)pSymmetric` | round-off | round-off | **1.9e-08 → 2.1e-13**, `O(h⁴)` | round-off |
| `GLRK(2)pSymmetric` | round-off | round-off | **1.9e-10 → 9.9e-15** | round-off |
| `GLRK(1)pMidpoint` | round-off | round-off | **1.3e-04 → 1.1e-07**, `O(h³)` | round-off |
| `GLRK(2)pInternal` | round-off | round-off | **2.1e-07 → 1.4e-12**, `O(h⁴)` | round-off |
| `GLRK(2)pLobattoIIIAIIIB` | round-off | round-off | **1.9e-10 → 6.9e-15** | round-off |

On `MasslessChargedParticle` **every one of the thirty Gauß-inner methods shows a clean
`O(hᵏ)` defect**, `k ∈ {3,4}`, with 100-step drifts up to 1.2e-05. Nothing is at round-off.
So the violated conditions do bite; they simply did not show on Lotka–Volterra.

#### Established: a linear `ϑ` makes any projection method symplectic for free

If **every** component of `ϑ` is at most linear in `q`, an unprojected variational
Runge–Kutta method already lands on the constraint manifold. With `ϑ(q) = Cq + d` and
`b̄ = b`,

```
ϕ(qₙ₊₁,pₙ₊₁) = h Σᵢ bᵢ [ (Cᵀ - C) Vₙ,ᵢ - ∇H(Qₙ,ᵢ) ] ,
```

which the stage equations annihilate. The projection is then inert, `Λ̃ = 0`, and every
`Λ̃`-carrying term in the wedge product drops out — so *any* projection built on such an
inner method is symplectic, whichever of the theorem's conditions it violates.

Verified directly (step 7): unprojected `VPRK` on the IODE form, `max|ϕ|` over 1000 steps
at `h = 0.1`:

| problem | `VPRKGauss(2)` | `VPRKGauss(3)` |
|:---|---:|---:|
| **`PointVorticesLinear`** — `ϑ = ½(-γ₁q₂, γ₁q₁, -γ₂q₄, γ₂q₃)`, linear | **2.6e-15** | **3.6e-15** |
| `PointVortices` | 2.2e-04 | 4.9e-05 |
| `LotkaVolterra2d` | 6.7e-02 | 2.3e-05 |
| `MasslessChargedParticle` | 1.7e-04 | 5.2e-07 |
| `MasslessChargedParticleSingular` | 3.1e-04 | 3.5e-07 |

And the corollary end-to-end (step 2): on `PointVorticesLinear` **all thirty** Gauß-inner
projection methods are at round-off — one-step defect `2e-16 … 4e-15` flat across
`h = 0.5 … 0.05`, 100-step drift `≤ 9e-15`, `max|ϕ| ≤ 4e-15` — including every construction
that violates one of the theorem's conditions. On the nonlinear `PointVortices` the same
thirty all fail badly (relative drift of order 1 at `h = 0.5`, `2.8e+00 … 4.2e+01` over 100
steps). The fitted exponents in that column are not meaningful there: the data is
non-monotone because the methods are diverging at the larger step sizes rather than
exhibiting a clean defect.

A test problem with a linear `ϑ` therefore certifies any projection method as symplectic
and carries no information about its tableau. This is the criterion to check before using
such results as evidence.

#### Open: why the other three problems also came out at round-off

An earlier revision of this section conjectured that an *affine second component* of `ϑ`
was the mechanism. That is a strictly weaker condition than the one above, and it is
**not** an instance of it — none of the four problems of the sweep is in the linear class:

| problem | `ϑ` | |
|:---|:---|:---|
| `LotkaVolterra2d` | `(q₂ + log(q₂)/q₁, q₁)` | `ϑ₁` nonlinear |
| `LotkaVolterra2dSingular` | `(log(q₂)/q₁, 0)` | `ϑ₁` nonlinear |
| `MasslessChargedParticleSingular` | `(-A₀x₂(1+2x₁²+⅔x₂²), 0)` | `ϑ₁` cubic |
| `MasslessChargedParticle` | `A₀/2(1+x₁²+x₂²)(-x₂, x₁)` | both cubic |

and step 5 shows the projection is genuinely active on all of them: `Λ̃` ranges from
2.4e-03 to 2.6e+01, never zero. So the round-off results on the first three remain an
unexplained empirical regularity, not a consequence of the linear-`ϑ` criterion. What is
established is the negative half: on `MasslessChargedParticle`, where no component of `ϑ`
is affine, every method drifts.

#### What the proof would need

At the point where the manuscript invokes `P̃ᵢ = ϑ(Q̃ᵢ)`, the `h¹` terms have been reduced to

```
Σᵢ b²ᵢ dG̃ᵢ ∧ dQ̃ᵢ + Σᵢ b⁴ᵢ dP̃ᵢ ∧ dΛ̃ᵢ
  = Σᵢ (b⁴ᵢ - b²ᵢ) dϑ(Q̃ᵢ) ∧ dΛ̃ᵢ  +  Σᵢ b⁴ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ ,     Φ̃ᵢ = ϕ(Q̃ᵢ,P̃ᵢ)
```

With `b⁴ = b²` the first sum dies and what is actually required is the **two-form**
condition

```
Σᵢ b⁴ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ = 0 .                                                (★)
```

The theorem discharges `(★)` by the far stronger pointwise demand `Φ̃ᵢ = 0`. Measuring both
factors (step 5) shows neither is the mechanism:

* `Λ̃ᵢ ≠ 0` for every method (1e-3 … 1e-8) — the projection is always active, so "the
  projection does nothing" is not the explanation.
* `Φ̃ᵢ = 0` to round-off for `pSymmetric` and `pInternal` at every `s`, but `Φ̃ᵢ ≠ 0` for
  `pSymplectic` (1.9e-04 … 1.8e-07) and for `pMidpoint` at `s = 2,3` — and *both* groups
  are at round-off on Lotka–Volterra.

For `pSymmetric` the pointwise condition is a **consequence, not an assumption**: its
projective stages are `(qₙ,pₙ)` and `(qₙ₊₁,pₙ₊₁)`, so its single `ω` row reads
`ϕ(qₙ,pₙ) + R∞ ϕ(qₙ₊₁,pₙ₊₁) = 0`, which propagates `ϕ = 0` by induction from a consistent
initial condition. The manuscript states this in its Symmetric-Projection example but does
not feed it back into the theorem. Restating the theorem with "the projective stages are
the endpoints and the initial condition is consistent" in place of `P̃ᵢ = ϑ(Q̃ᵢ)` would make
that hypothesis checkable from the tableau rather than from the solution.

That still leaves `pSymmetric`'s `cond4` violation. With `σ = 2`, `ρ = 1` the `δ` row ties
`Λ̃₁ = R∞Λ̃₂` and the `ω` row ties `Φ̃₁ = -R∞Φ̃₂`, so the multiplier block collapses to one
dimension and `(★)` reduces to `R∞(b⁴₁ - b⁴₂) dΦ̃₁ ∧ dΛ̃₂`, which vanishes for `b⁴₁ = b⁴₂`,
i.e. `R∞ = +1`. *(Derived here, not verified numerically.)* A rank argument of this kind —
the `δ` constraints collapsing the multiplier space — looks like the most promising way to
weaken the hypotheses, but it does not by itself explain the `MasslessChargedParticle`
results, where all of these methods do drift.

### Option 3 — composed symplectic projection

`TableauVSPARKSymplecticProjection` *is* the composed construction: a symplectic Gauß
inner method plus a projection carried by separate Lobatto projective stages, with the
endpoint constraint in the last row of `ω_λ`. It escapes the S16 no-go because it does
**not** impose `ϕ = 0` at every internal stage — only on the projective stages — so the
counting argument that forces `μ` never arises, and the constraint is still held at
`1e-15`. But per S17 it is *not* symplectic on a one-form with no affine component, and it
is the least robust member of the family under a change of gauge, so it does not resolve
the question that S16 leaves open.

What remains genuinely unexplored: the VSPARK-primary family projects on the **primary**
constraint only. SLRK's distinguishing feature — averaging the **secondary** constraint —
has no composed-projection counterpart yet. A "VSPARK-secondary with symplectic
projection" is the open construction; `VSPARKsecondary` exists but uses the
constraint-at-every-stage layout that S16 rules out.

# Release Notes

All notable changes to GeometricIntegrators.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually changed,
so that a compat-only bump can be told apart from a rename or a change in results.

Entries for 0.15.0 and later were written from the git history; the 0.4.0 – 0.11.0 entries below
are the original release notes, kept verbatim. Versions 0.12 – 0.14 were never written up and
remain a gap.


## 0.18.3

### New Features

* **`CGVINodal`, a continuous Galerkin variational integrator on a nodal basis.** Where
  [`CGVI`](@ref) imposes the continuity of the trajectory across the interval boundaries
  weakly, through Lagrange multipliers added to the action, `CGVINodal` builds it into the
  basis, following [OberBloebaum:2015](@cite). The basis must be interpolatory with nodes at
  both ends of the interval — a Lagrange basis on Lobatto-Legendre nodes — so the
  coefficients *are* nodal values: the first is pinned to the known `qₙ`, the last *is*
  `qₙ₊₁`, and the momentum is computed explicitly rather than solved for. That leaves
  `D*(S-1)` unknowns against `CGVI`'s `D*(S+1)`.

  Measured convergence order is 2s-2 in the position on Lobatto(s) nodes, the same as
  `CGVI` on Gauss(s), and the two are of comparable accuracy at equal `s`.

  The constructor rejects a basis whose nodes do not include the endpoints. The formulation
  reads `q(0)` off `X[1]` and `q(1)` off `X[S]` rather than reconstructing them, so a Gauss
  basis builds fine and then integrates a wrong trajectory — the same class of silent error
  `DGVIP0` already asserts against.

  The integrator comes from NonlinearIntegrators.jl, where it was the linear reference the
  network integrators were compared against, and is removed there in favour of this one.

### Changes

* **`CGVI` and `CGVINodal` share one core**, `src/integrators/cgvi/integrators_cgvi_common.jl`,
  in the shape the DGVI family already uses: an abstract `CGVIMethod` supertype carrying the
  traits, the solver defaults, the accessors, `cgvi_coefficients`, `show`, the cache, the
  initial-guess helpers, `components_q!`/`components_v!`/`components_p!`, and the
  `residual!`/`update!`/`integrate_step!` entry points. Each variant supplies only
  `description`, `solversize`, `initial_guess!`, `components!`, `residual!` and `update!`.

  `CGVI`'s results are unchanged, bit for bit.

* `CGVICache` is now shared by both variants. It lost the `s̃` field, which was allocated and
  never read, and takes the number of degrees of freedom and the solver size as constructor
  arguments. It is not exported and no user code should name it. `CacheType` deliberately
  stays a function of the method alone: `CacheDict` type-asserts its `getindex` on it, and a
  `CacheType` that reads a value off the problem does not constant-fold, which leaves the
  cache inferred abstractly and makes the Newton hot path box on every stage access.

* Both variants now define `solversize`, as every other implicit method here does. That is
  what sizes the nonlinear solution vector and what `initsolver` scales the default
  `f_abstol` by; previously the CGVI family took the `solversize == 0` fallback.

### Documentation

* `docs/src/integrators/cgvi.md` gains a section on the nodal formulation and the
  [OberBloebaum:2015](@cite) reference.


## 0.18.2

### Breaking Changes

* Requires SimpleSolvers 0.12.1 and GeometricIntegratorsBase 0.6.3. Coupled as before:
  GeometricIntegratorsBase 0.6.3 requires SimpleSolvers 0.12.1, so neither range is satisfiable
  with the old pin on the other.

  Of the three breaks in SimpleSolvers 0.12, one reaches this package: **a `NonlinearSolver` no
  longer emits line-search warnings from inside its iteration**. `solver_step!` stopped calling
  `linesearch_warnings`, so a rejected line search now reports to its *caller* through the returned
  status and to nobody at all if the caller drops it. Every solve here dropped it, because every
  solve here called `solve!`. That is what the change below is for. The other two — a new field on
  `NonlinearSolverStatus`, and a `LinesearchMethod` implementing `solve_with_status` rather than
  `solve` — do not reach a package that constructs neither.

  `Bisection`'s two fixes in 0.12 (it no longer claims success when it cannot bracket, and it
  bisects toward a minimum rather than toward whichever root the value bracket held) reach no line
  search here: `default_linesearch` is `Backtracking` and every method solves with `Newton()`. The
  one place this package does reach for `Bisection` is not a line search at all — PGLRK bisects the
  *energy residual* over λ through `SimpleSolvers.bisection`. That call **does** take the
  no-bracket branch, on every step of every `λmax = 1E-300` run in the test suite, and it is
  unaffected all the same: the endpoint returned there is bit-for-bit the one 0.11 returned
  (`abs(y₀) ≤ abs(y₁) ? α₀ : α₁`), the warning that accompanies it is gated behind the
  `verbosity = 0` that `solve_λ!` passes and predates 0.12 anyway, and the success flag that 0.12
  corrected is the one thing `bisection` does not hand back to its caller — `solve_λ!` re-derives
  the failure itself, from `abs(λ) ≥ λmax`.

### New Features

* **Every nonlinear solve goes through `solve_with_status!`**, and the status it returns is handed
  to GeometricIntegratorsBase 0.6.3's `check_solver_status`. Nineteen call sites: fourteen
  `integrate_step!` methods across the RK, VI, CGVI, DGVI, DVI, HPI and SPARK families, PGLRK's
  `energy_residual!` — a stage solve *inside* a step rather than a step, which is why it is
  counted apart — the three projection integrators, and DIRK's per-stage loop. It replaces the
  pair of commented-out stubs
  (`# println(status(solver))` / `# println(meets_stopping_criteria(status(solver)))`, and the
  older `# print_solver_status(int.solver.status, …)` spelling) that had sat under the solve in
  most of those files.

  The fifteen that own a persistent `solverstate(int)` use the state-taking form of
  `solve_with_status!`, added in SimpleSolvers 0.12.1 at GeometricIntegratorsBase's request; the
  other four have no state to reuse and take the state-building form. DIRK's per-stage solvers get
  a `NullSolverState` from the `SingleStageSolvers` wrapper, and a `ProjectionIntegrator` has no
  `solverstate` field at all.

  `check_solver_status` is silent by default — SimpleSolvers remains the one voice that reports a
  failed solve, and warning here as well would say the same thing twice per time step. So **no run
  changes what it prints**. See GeometricIntegratorsBase's notes for what the hook is for and how
  to override it, and for why reading the status now costs nothing per solve.

* PGLRK is the one method where the status needed a decision rather than a rewrite. Its stage
  solves happen inside `energy_residual!`, which `SimpleSolvers.bisection` calls once per trial λ
  — including at the endpoints `±λmax`, where a solve that struggles is expected rather than
  exceptional. Passing every probe's status to `check_solver_status` would hand a caller who had
  overridden it a bisection probe to reject. `solve_λ!` checks the *accepted* λ once instead,
  reading its outcome off the persistent solver state that the final `energy_residual!` call
  leaves behind — which costs nothing and cannot describe a different solve than the one the step
  keeps. This is the one place here that still calls `status(s, state)` by hand, and it is sound
  because SimpleSolvers 0.12.1 states the property it relies on outright: nothing touches the state
  between the end of the solve loop and the caller, so rebuilding the status from it afterwards
  gives the same value. Its own test suite pins that.

* The seven `solve!` calls under `src/integrators/vprk/` are deliberately left alone. That whole
  directory is dead: it is reached only through `src/integrators/VPRK.jl`, whose every `include` is
  commented out and which is itself never included. Converting it would imply it is live.

### Tests

* `test/verification/pglrk_convergence_tests.jl` counts the statuses PGLRK hands to
  `check_solver_status` and requires exactly one per time step. This is the one decision in this
  release that no other assertion could catch: the hook's default returns its argument, so calling
  it once per *bisection probe* instead of once per step is invisible to every other test in the
  suite, and would only surface for a caller who had overridden it to reject a non-converged step —
  who would then find bisection probes being rejected.

  `solve_λ!` has **three** ways through it, each with its own `check_solver_status` call site, and
  the test drives all three at 5 checks over 5 steps. Which route a run took is *asserted* and not
  assumed, from the count of fallback `@debug` lines captured off a `TestLogger`, because a hook
  count of `nsteps` is what every route is supposed to produce and so cannot by itself tell them
  apart:

  | route | problem, `λmax` | fallbacks |
  |:---|:---|:---|
  | bisection locates a root | Lotka-Volterra, default | 0 |
  | bisection finds no sign change, step falls back to plain Gauss | Lotka-Volterra, `1E-300` | 5 |
  | early return, the unperturbed method already conserves the energy | harmonic oscillator, `1E-300` | 0 |

  The third row needs a problem where plain Gauss conserves `h` to the `ftol` of the early return,
  which the harmonic oscillator is and Lotka-Volterra is not: Gauss preserves the quadratic
  invariant of a linear system exactly (|h − h₀| ≈ 7E-18 against ftol ≈ 1.8E-15), whereas on
  Lotka-Volterra at Δt = 0.1 it drifts to ≈ 6E-11 against ftol ≈ 3.6E-15, four orders too far. Note
  that `λmax` does *not* reach the early return, whose condition is on the λ = 0 residual alone —
  which is why the third row keeps the `λmax` of the second: had the early return not been taken,
  the bisection would have run in [±1E-300] and fallen back on every step exactly as it does in row
  two, so `nfallback == 0` there is what says the bisection never ran at all.

* **`VSPARK(SPARKLobattoIIIBIIIA(2))` accepts either outcome, which is what CI on Julia 1.13 and
  nightly was failing on.** That case raises `SingularException: Zero pivot found at index 25` on
  1.13 and nightly under Linux and Windows, while returning a ~1E-6 answer on 1.10 and 1.12
  everywhere and on macOS throughout. The failure is older than this release — it is on `main` at
  the previous commit with the identical signature — and it is not a regression to chase: the stage
  system is **numerically singular**, cond ≈ 5.6E16 with σmin = 5.7E-17 against σmax = 3.2, which is
  the same matrix to within a factor of 1.5 as the `SPARKLobattoIIIAIIIB(2)` and `SPARKGLRK(2)`
  siblings the suite already asserts `SingularException` for. Its σmin sits an order below the
  `n·eps·σmax ≈ 1.8E-14` at which a 26×26 system stops having a numerical rank.

  Whether LAPACK's `getrf` lands on an exact zero pivot or on one of ~1E-17 is a rounding accident
  of the BLAS kernel, so the outcome is genuinely platform-dependent. The `@test` on this case was
  promoted from `@test_broken` under SimpleSolvers 0.10 on the strength of one platform where it
  survived; that promotion is retracted, and the test now asserts the answer where there is one and
  `err isa SingularException` where there is not. Audit finding S8 always covered the case — see the
  S8 retraction in `docs/src/audit.md`, which carries the singular-value measurements.

### Documentation

* **`docs/src/audit.md` carries the re-measured warning census.** This closes the item that stood
  under `## Open Issues` as *"`docs/src/audit.md` still carries the pre-0.12 warning census"*: the
  page now has its own `Update (SimpleSolvers 0.12.1 / GeometricIntegratorsBase 0.6.3)` block with
  the 14 → 8 figures recorded under *Findings* below, and with the one thing 0.12 changes about the
  reasoning already on that page — the `maxlog` caps having become a back-off, the arithmetic by
  which three saturating sites gave exactly 9 warnings no longer describes the mechanism, and
  neither does the `@test_nowarn`-goes-blind caveat that followed from it.

* The same page records the **S8 retraction** for `VSPARK(SPARKLobattoIIIBIIIA(2))` described under
  *Tests* above, with the singular values behind it, and strikes the 0.17.0 promotion where it was
  originally recorded rather than leaving two entries that disagree.

### Findings

Test outcomes are unchanged across this bump, and this was measured rather than assumed: the suite
was run at `HEAD~` against SimpleSolvers 0.11.0 and GeometricIntegratorsBase 0.6.2 in a separate
worktree, and all **28 testsets report identical pass and broken counts** — the six route
assertions added above excepted, since they do not exist there.

**The warning census moved, and this is what it moved to.** A full run emitted 14 warnings before
and emits 8 now. Unchanged: the seven `Nonlinear solver failed at timestep n=…: non-finite
direction vector` (GeometricIntegratorsBase's own, at n=5, 7, 9×3, 10×2) and the one stagnation
warning from `VSPARK(SPARKLobattoIIIBIIIA(2))` — the one entry here that is platform-dependent,
since it is emitted exactly on the runs where that solve returns rather than raising
`SingularException`, as the *Tests* entry above describes. Gone: **six `Backtracking line search:`
warnings** —
three `no step satisfied the sufficient decrease condition in 13/14 trials` and three
`φ'(0) = <tiny positive> (with φ(0) = …)`. Those are exactly the per-iteration line-search warnings
SimpleSolvers 0.12 stopped emitting from inside `solver_step!`.

**Nothing was lost with them.** The stagnation warning reports the same residual to the bit
(`rfₐ = 3.571488053775618e-5`) and the same diagnosis, and 0.12 appends the cause to it:

> The line search reported `LINESEARCH_NO_DESCENT` on 2 of the 3 step(s), i.e. φ'(0) > 0 — the
> direction was not a descent direction at all, which points at the Jacobian rather than at the
> tolerance: a stale one under `refactorize > 1`, a nonzero `regularization_factor`, or an inexact
> linear solve.

That is the same three `φ'(0) > 0` events the deleted per-iteration warnings were reporting one at
a time, now counted and named once at the end of the solve. `docs/src/audit.md`, where the census
lives, carries these numbers as its own update — together with the one thing they change about the
reasoning there: SimpleSolvers 0.12 replaced the `maxlog` caps with a back-off (occurrences 1, 2,
4, 8, …), so the "three `maxlog = 3` sites saturate at exactly 9" arithmetic of the earlier census
no longer describes the mechanism. This suite never reaches the back-off regime, the repeating
diagnosis now firing once per run, so it is inherited but not exercised from here.


## 0.18.1

### Breaking Changes

* Requires SimpleSolvers 0.11 and GeometricIntegratorsBase 0.6. The two bumps are coupled:
  GeometricIntegratorsBase 0.6 requires SimpleSolvers 0.11, so neither range is satisfiable with
  the old pin on the other. Nothing this package calls was removed by either.
* **A non-finite direction now aborts the integration.** SimpleSolvers 0.11 tests
  `all(isfinite, …)` where 0.10 tested only `isnan`, so an *overflowed* direction — which used to
  pass every guard and stall silently, `Inf * nan_factor` being `Inf` — raises
  `NonlinearSolverException`. GeometricIntegratorsBase catches it in the time-stepping loop,
  warns naming the timestep, and breaks before the `copy!`, so the solution holds valid data up
  to `n-1` and zeros after it.

  That warning is GeometricIntegratorsBase's own: the `verbosity = 0, warn_iterations = 0` pair
  the SPARK suites pass does not reach it, and it carries no `maxlog`, so its count is the true
  number of events. Seven fire per test run, all from methods already recorded `@test_broken` for
  diverging. They are left unsuppressed deliberately — each is a correct symptom, and a more
  informative one than the silent stall it replaces, since it names the step at which the
  trajectory stops meaning anything.

* `default_options` inherits `f_stall_window = 50` from GeometricIntegratorsBase, SimpleSolvers
  0.11's criterion for a solve whose residual sits on a floor above `f_abstol` while the iterate
  keeps moving — the case `max_stalls` cannot see. No measured effect on this test suite. Options
  are merged rather than replaced, so a method needing a longer window overrides it with one
  keyword.

### New Features

* DIRK calls the typed `default_linesearch(eltype(x), method)` that GeometricIntegratorsBase 0.6
  added alongside the untyped hook.

Test outcomes are unchanged across this bump: all 24 testsets report identical pass and broken
counts, and Aqua passes.


## 0.18.0

### Breaking Changes

* **Four Runge-Kutta method types are renamed**, to avoid a collision with the methods of the
  same name that GeometricIntegratorsBase 0.5.2 now exports. `Integrators.jl` reexports that
  package, and two distinct types sharing a name resolve to nothing once both are in scope. The
  Runge-Kutta side takes the suffix, following `ImplicitEulerRK`, which was renamed the same way
  and for the same reason:

  | old | new |
  |---|---|
  | `CrankNicolson` | `CrankNicolsonRK` |
  | `ImplicitMidpoint` | `ImplicitMidpointRK` |
  | `SymplecticEulerA` | `SymplecticEulerARK` |
  | `SymplecticEulerB` | `SymplecticEulerBRK` |

  Only the Julia types are renamed. The tableaus keep their names, so `TableauCrankNicolson` and
  `TableauImplicitMidpoint` are untouched, the partitioned tableaus built inline stay
  `:SymplecticEulerA` and `:SymplecticEulerB`, and the bibliography keys behind
  `reference(Val(:CrankNicolson))` and `reference(Val(:ImplicitMidpoint))` are unchanged — those
  name the methods, not the wrappers. The unsuffixed names now resolve to
  GeometricIntegratorsBase.

* **CGVI and DGVI results change.** CompactBasisFunctions 0.3 fixes the Legendre derivative,
  which had been missing the √(2j+1) normalisation of its own basis functions.
* Requires RungeKutta 0.6, which dropped the `get_` prefix from its accessors and then removed
  the node and weight ones altogether, having become one-line forwardings to QuadratureRules.
  The quadrature data is now taken from QuadratureRules directly, which is one indirection fewer
  and needs no new dependency. `lobatto_nullvector` stays and is imported explicitly, in both
  `Integrators` and `SPARK`, as it is no longer exported.
* Requires QuadratureRules 0.2, CompactBasisFunctions 0.3, GeometricBase 0.14.8 and
  GeometricIntegratorsBase 0.5.2.
* `basis` and `nnodes` are extended from GeometricBase rather than defined here. `basis(::CGVI)`,
  `basis(::DGVIMethod)` and their `nnodes` counterparts were bare definitions, so they created
  GeometricIntegrators-local functions rather than extending the ones CompactBasisFunctions and
  QuadratureRules use — and distinct functions sharing a name resolve to nothing when two of them
  are in scope together. `nbasis` stays imported from CompactBasisFunctions, which owns it.
* The `GenericLinearAlgebra` dependency is dropped. It was never called by name; the `using` in
  `src/SPARK.jl` was left over from 2020, when tableau construction still lived here and needed
  generic `eigvals` for BigFloat companion matrices. One call did rely on it, in the DVI
  integrators rather than in SPARK: `rank(Matrix(tableau.a))` routes through `svdvals`, and there
  is no stdlib `svdvals!(::Matrix{BigFloat})`, so `DVRK(TableauGauss(BigFloat, s))` only worked
  because RungeKutta happened to load GenericLinearAlgebra for its own use. The check is a
  diagnostic `@warn`, so the rank is now taken in Float64.

  Motivation: GenericLinearAlgebra 0.4.0 cannot precompile on Julia 1.13, where overwriting its
  `LinearAlgebra.eigencopy_oftype` method for `UpperHessenberg` is a hard error.

* Documenter, SharedArrays, Test and Parameters are removed as dependencies; none is used by the
  loaded sources. The inert `[extras]` block goes with them — there is no `[targets]` section, so
  it has had no effect since the move to a workspace layout. `ProgressMeter` is deliberately
  kept, as a marker for the commented-out simulation code.

### Fixes

* **`R∞` operator precedence in sixteen Lobatto-pair projection constructors.** They passed
  `R∞=-1^(s+1)`, which Julia parses as `-(1^(s+1))` and is therefore −1 at every s, where the
  intended `(-1)^(s+1)` is +1 at odd s. The Gauss siblings in the same functions already wrote
  `(-1)^s` with parentheses, so the two forms sat side by side in one file. `R∞` sets β,
  d_λ, δ_λ and β_q/β_p, and at odd s the wrong sign zeroes β. All sixteen change at odd s, none
  at even s (even-s results are bit-identical). Two methods change order:

  | method | order | error |
  |---|---|---|
  | `VSPARKLobattoIIIBIIIApSymmetric(3)` | measurement failed → 4.000 | 4.39e-07 → 1.89e-11 |
  | `VSPARKLobattoIIIAIIIBpSymmetric(3)` | 3.094 → 3.999 | 4.17e-11 → 1.21e-11 |

  `VSPARKLobattoIIIBIIIApSymmetric(3)` is promoted from broken to a real convergence assertion.
  The finding previously filed against it in the verification report as inherent divergence is
  retracted and re-filed as a fixed bug.

* The Vandermonde solve of `lobatto_gauss_coefficients` is pinned to BigFloat. The RungeKutta 0.6
  migration had to spell this out once the quadrature data started coming from QuadratureRules,
  whose accessors default to Float64 rather than BigFloat; solving in Float64 instead moves the
  coefficients by 2.0e-15 in relative norm at s = 4 and by 9.3e-13 at s = 8, and no test would
  have caught it, since the result is converted to Float64 either way.

### Tests

* Aqua's stale-dependency check is enabled, now that the dead `[compat]` entries are gone. It
  would not have caught Documenter or ProgressMeter, both of which have a live `using`.


## 0.17.0

### Breaking Changes

* Drops GeometricIntegratorsBase 0.3/0.4 and SimpleSolvers 0.9 from compat, so environments that
  resolved against those can no longer resolve GeometricIntegrators.
* **The framework `f_abstol` default changes** from a flat `8eps()` to
  `max(8, solversize(method, problem)) * eps(datatype(problem))`, which comes from
  GeometricIntegratorsBase 0.5.1. This shifts the machine-precision error floors of the
  high-dimensional SPARK methods, whose `f_abstol` now resolves to 1.8e-15 … 1.1e-14 rather than
  1.8e-15 throughout. The affected test thresholds are re-measured.
* `solversize` and `nullvectorsize` take their arguments in `(method, problem)` order, following
  GeometricIntegratorsBase 0.5.1 — 21 definitions plus every call site. The two ad-hoc `f_abstol`
  overrides this package carried (a flat `8eps()` in `src/Integrators.jl` and SPARK's `8e-15`)
  are removed in favour of the framework default.

### Fixes

* `HSPARKsecondary` warning suppression restored on the `@test_throws` loop, where it had been
  dropped rather than substituted when those cases were rewritten from `@test_broken`.

### Tests

* The `muffle` log-suppression scaffolding is removed throughout — six definitions in two
  variants across 42 call sites — in favour of `verbosity = 0, warn_iterations = 0` on the
  genuinely divergent cases. SimpleSolvers 0.10 makes this possible: a line search now shares its
  solver's `Options`, so `verbosity` finally reaches it. `verbosity = 0` silences only
  SimpleSolvers, where a `NullLogger` also swallowed unrelated failures.
* `@test_nowarn` tripwires are added around the converging variational and Hamilton-Pontryagin
  convergence calls, now that which orbit warns no longer depends on platform floating-point
  details.
* The SPARK `@test_broken` cases now assert the failure mechanism where it is unambiguous:
  structurally singular stage systems assert `@test_throws SingularException`, order-deficient
  and unstable methods assert `converged`, stalling solves assert `stalled`. Marginally singular
  cases, which converge or zero-pivot depending on rounding, stay `@test_broken`.
* `VSPARK(SPARKLobattoIIIBIIIA(2))`, formerly a singular `@test_broken`, converges under
  SimpleSolvers 0.10 and is promoted to `@test`. See Open Issues — its solve still stalls well
  above the tolerance it was asked for.
* Every magnitude and error tolerance is normalised to the tightest {1,2,4,8}×10ⁿ strictly above
  the freshly measured value; dropping the "5" mantissa makes the powers-of-two grid guarantee
  measured < TOL ≤ 2 × measured.


## 0.16.10

### Breaking Changes

* Requires GeometricProblems 0.8. Every dependency gains a compat entry, and the test and docs
  dependencies are bounded where they are resolved.

### Fixes

* **SLRK null-vector term.** The null-vector multiplier was added to both the momentum-stage row
  and the primary-constraint row. The first lives in Z-space (P = p + h·Z) so it carries an extra
  factor h, and the two combined to (1-h)·μ·d_i/b_i. The stage Jacobian was therefore exactly
  singular at Δt = 1 and ill-conditioned near it — cond(J) 1.5e2 → 2.2e3 → 2.3e4 → 1.6e11 as
  Δt → 0.9, 0.99, 0.999, 1.0, now flat at ~47. No other SPARK integrator does this. Present since
  the original implementation, and behaviour-neutral at usable step sizes: trajectories unchanged
  to ≤ 2 ulp.
* `SLRKLobattoIIICbarC` was registered under the same name symbol as `SLRKLobattoIIICCbar`.
* The SLRK docstring was a verbatim copy of the `VSPARKsecondary` one and described a different
  scheme. Rewritten from the code and the manuscript, including the requirement that the
  `LDAEProblem`'s `f` be the full dL/dq — the opposite convention to `VSPARKsecondary`, previously
  undocumented and a silent source of wrong results.
* **DVRK traits.** `order`, `issymmetric` and `issymplectic` called an undefined
  `tableaus(method)` and threw `UndefVarError` on every `DVRK` instance. This is why the method
  table reported `missing` for DVRK and why its convergence tests had to pass `expected =`
  explicitly.
* **The DVRK tableau invertibility check false-positived on well-conditioned tableaus.** `det` is
  not scale-invariant and the determinant of an RK coefficient matrix shrinks rapidly with the
  number of stages, so `abs(det(a)) ≤ eps^(3/4)` flagged `Gauss(10)` and up as singular
  (det = 1.5e-12 at cond = 1.2e2). Replaced by `rank(Matrix(a)) < s`, which separates Gauss (full
  rank) from LobattoIIIA/IIIB (rank s-1, cond ≥ 6e17) cleanly and works for Float32 and BigFloat
  tableaus too.
* **The even-dimension guard is widened.** `DVIA`, `DVIB`, `CMDVI` and `CTDVI` fail on odd D
  exactly as `DVRK` did — `div(D,2)` dropped a component silently, and for D = 1 left the
  nonlinear system underdetermined. The guard moves to `check_dvi_dimension` and is called from
  all three cache constructors, so all five methods now give the same `ArgumentError` instead of
  a `SingularException`.
* Minor fix in the initial guess for `CTDVI`; the SPARK `show` methods are fixed and the show
  tests now assert their titles.
* `solversize` becomes one function with one definition.

### Documentation

* **DVRK's order-2s claim holds only within a gauge class.** DVRK requires L = ϑ(q)·q̇ − H(q)
  with d/2 components of ϑ vanishing identically; ϑ is only defined up to an exact one-form, and
  the gauge decides whether a given problem satisfies that hypothesis. Every existing DVRK test
  used `LotkaVolterra2d`, whose gauge has ϑ₂ ≠ 0 and is therefore outside the class. In class
  `DVRK(Gauss(s))` attains the full order 2s; out of class it stays convergent at s. The same
  split shows up in symplecticity: the defect sits at the finite-difference floor in class, and
  grows with h out of class. The "order reduction on degenerate Lagrangians" previously recorded
  as a property of the method was a testing artefact.
* The hypotheses under which DVRK is symplectic are documented and checked: the tableau
  symplecticity condition and an invertible coefficient matrix, an even-dimensional configuration
  space, and p₀ = ϑ(q₀). `DVRK(RadauIIA(s))` and `DVRK(LobattoIIIA(s))` now warn rather than
  silently producing a non-symplectic method. Note that the vanishing-components condition and
  the invertibility of ∂ϑ_μ/∂q^ν are *separate* conditions; earlier prose folded one into the
  other.
* A full SLRK section in `docs/src/integrators/spark.md`: scheme, ω construction, null vector,
  constructor table, gauge caveat, and what the family does and does not preserve.
* The verification report is published as a documentation page (`docs/src/audit.md`) rather than
  living in the repository root.
* `scripts/` gains `dvrk_convergence.jl`, `dvrk_symplecticity.jl`,
  `projected_vprk_symplecticity.jl`, `slrk_verification.jl` and
  `vspark_projection_symplecticity.jl`, which reproduce every number in the report.

### Findings

Three negative results, recorded in `docs/src/audit.md` and carried into Open Issues below:

* **SLRK is not symplectic.**
* **The VSPARK projection methods are not empirically symplectic either.**
* **The midpoint projection is not symplectic**, and `VPRKpSymplectic` is the same map as
  `VPRKpStandard`.

### Continuous Integration

* Coverage is computed in a single job rather than in all twelve matrix jobs.


## 0.16.9

### New Features

* `IRK3`, a fully implicit Runge-Kutta method.


## 0.16.8

### Breaking Changes

* Requires GeometricProblems 0.7.
* The five DGVI variants share one abstract method type, one superset coefficient block, one
  pruned cache and a common integrator layer, which removes 1714 lines net while taking the
  family from zero working variants to five.

### New Features

* DGVI, FLRK and PGLRK are ported to the method-based `GeometricIntegrator` architecture,
  registered as methods, documented, and covered by the same static-plus-dynamic verification
  standard as the rest of the package. `Discontinuity` and `PathIntegral` are documented.

### Fixes

Nine defects, four of them numerical:

* `CoefficientsPGLRK` could not construct at all: `new(...)` passed 11 arguments to a 14-field
  struct, since `@CoefficientsRK` gained â/b̂/ĉ, so P landed in â and the Q matrix in
  `b̂::Vector`. `Legendre` was also never imported.
* FLRK took its stage force from the projection field g = (∇ϑ)ᵀv instead of f = (∇ϑ)ᵀv − ∇H, so
  the momentum update was missing −∇H. Reverting the fix degrades max|p − ϑ(q)| from 1E-15 to
  9.9E-2 while leaving q correct, which is why the old q-only test never caught it.
* FLRK as an `RKMethod` would have been silently replaced by a plain IRK by
  `initmethod(::RKMethod, …)`; it is now an `LODEMethod`.
* `VPRKpTableau` referenced unbound names and needed a non-dependency (NLsolve); rewritten with
  the multipliers folded into one coupled Newton system.
* `CoefficientsPGLRK(2)` is inconsistent by construction and is now rejected.
* `VPRKpTableau` needs s ≥ D+1 to fit its multipliers and s ≥ D+2 for full order.
* `DGVI`/`DGVIP0`/`DGVIP1` are capped at order 2⌊s/2⌋ by their nodal-value flux.
* The disabled DGVI tests used a regular Lagrangian, on which DGVI's closure row collapses to
  q − p and the Jacobian is singular; they now use degenerate ones.
* The PGLRK symplecticity claim is corrected.

Verified afterwards: FLRK reaches order 2s in both q and p; PGLRK conserves energy to 3.6E-15
against Gauss's 7.3E-10 on a nonlinear problem; `VPRKpTableau` enforces the Dirac constraint to
1E-15 and recovers the order plain VPRK loses on degenerate Lagrangians; DGVIPI and DGVIEXP reach
the full order 2s.


## 0.16.7

### Fixes

* **`McLachlan4` was order 2, not 4.** The coefficients were composed as φ then adjoint, where
  McLachlan 1995 eq. 1.8 requires the adjoint first. The coefficients themselves are correct.
  The existing accuracy tolerance tightens 5E-4 → 5E-8.
* **HSPARK-secondary momentum coefficients** `a_p_2`/`a_p_3` were s×s but indexed over the σ
  projective stages, giving a `BoundsError`. Rebuilt as the s×σ conjugate-symplectic partners of
  `α_q_2`/`α_q_3`.
* **HSPARK-secondary null-vector** residual and component code was commented out while the cache
  still allocated the μ unknown, leaving an unconstrained Jacobian row and a `SingularException`.
  Re-enabled, matching the working VSPARK-secondary. A deeper singularity remains in its ω
  secondary-constraint block, so the EXPERIMENTAL family stays `@test_broken`.
* **HSPARKsecondary initial guess** built its NamedTuple with `v`/`f` rather than the solution
  derivative fields `q̇`/`ṗ` that every other SPARK integrator uses, so `solutionstep!` could not
  consume it (`FieldError` on hdae problems).
* δ-constraint residual row index R-1 → i in `integrators_hspark.jl` and
  `integrators_hspark_primary.jl` — neutral for P = 1, correct for P > 1.
* The momentum (p) error checks were NaN because `relative_maximum_error(…).p` divided 0/0 at the
  harmonic oscillator's p(0) = 0. Fixed upstream in GeometricSolutions 0.6.5; the checks are
  re-enabled.
* `KraaijevangerSpijker` is order 1, not 2 (RungeKutta 0.5.22 corrects the attribute; the method
  genuinely satisfies only the order-1 conditions, Σᵢbᵢcᵢ = 2 ≠ 1/2).
* The VPRK docstring's symplecticity condition is corrected to the general barred form.

### New Features

* `Yoshida6` (6th order, SS m=7) and `Yoshida8` (8th order, SS m=15) symmetric composition
  methods, extending the splitting methods beyond 4th order. Coefficients verified to high
  precision against Yoshida 1990 and McLachlan 1995 Table 2.

### Documentation

* New pages for the degenerate variational integrators (`dvi.md`: DVIA/DVIB/CMDVI/CTDVI + DVRK)
  and the Hamilton-Pontryagin-Galerkin framework (`hpg.md`), the latter re-enabled in `make.jl`;
  the CGVI page is enabled too. `spark.md` is expanded with the construction and the
  order/stability caveats. Stage-time and summation-index typos fixed in `rk.md`.

### Tests

* `test/verification/`: a reusable convergence and energy harness plus per-family convergence
  tests (rk, prk, splitting, variational, galerkin, dvi, hpi, projection, spark), wired into
  `runtests.jl` — 64 SPARK order assertions and 11 recorded deficiencies among them.
* Roughly 30 previously-commented VSPARK/HSPARK cases are re-enabled with honest tolerances, and
  the genuinely broken ones recorded as `@test_broken` with root causes rather than tuned away.
  Most SPARK non-convergence is an inherent method property, not a bug: symplectic plus
  constraint-at-solution reduces order or diverges; R(∞) = (-1)^(s+1) ≠ 1 drops GLVPRK and
  HPARKGLRK from 2s to 2 at s = 2; coinciding tableau pairs give a singular stage system at s = 2.
* Error bounds more than an order of magnitude looser than the measured error are tightened, from
  measurement.


## 0.16.6

### Breaking Changes

* Updates for GeometricIntegratorsBase 0.4.0, whose `SolutionStep` `reset!` takes a time rather
  than a timestep.


## 0.16.5

* GeometricProblems 0.7 in the test extras.


## 0.16.4

### Breaking Changes

* Requires SimpleSolvers 0.9: `NewtonMethod` is replaced by `Newton`, following the new naming
  upstream, and the default linesearch becomes `StrongWolfe`.

### Tests

* Some problematic SPARK tests are disabled and the remainder switched to the default solver
  settings.


## 0.16.3

### Breaking Changes

* Adapts to the removal of the dimension type parameter in `Cache`
  (GeometricIntegratorsBase 0.3); requires Parameters 0.12 or 0.13.

### Fixes

* An `IntegratorCacheSPARK` type parameter mismatch in `integrators_spark_parameters.jl`.
* The `cache.t̄`/`q̄`/`p̄` assignments are removed from `GeometricBase.reset!` for
  `IntegratorCacheSPARK`.

### Tests

* Tests for the Radau and Lobatto integrators.


## 0.16.2

### Documentation

* Cleanup in the documentation.


## 0.16.1

### Documentation

* Documentation fixes.


## 0.16.0

### Breaking Changes

* Requires GeometricBase 0.14, GeometricEquations 0.21, GeometricSolutions 0.6,
  GeometricIntegratorsBase 0.2 and PrettyTables 2 or 3.
* Runge-Kutta methods account for the precision of the solution, and the tableaus become
  type-dependent.

### Documentation

* Plots in the documentation are replaced by Makie.

### Tests

* Test dependencies are reorganised into workspaces, and `test/Project.toml` is checked in.


## 0.15.5

### Breaking Changes

* Requires SimpleSolvers 0.7.5 and the corresponding GeometricIntegratorsBase release.


## 0.15.4

### Breaking Changes

* Admits PrettyTables 3 alongside 2.

### Continuous Integration

* CI updated for Julia 1.13.


## 0.15.3

### Breaking Changes

* Requires SimpleSolvers 0.7.


## 0.15.2

### Breaking Changes

* Requires SimpleSolvers 0.6.

### Fixes

* Imports and dependencies cleaned up; a compat entry added for GeometricProblems; the missing
  `base.md` added to the documentation.


## 0.15.1

### Fixes

* Adapts to the initial guess fixes in GeometricIntegratorsBase, and fixes the initial guess for
  the IRK implicit method.
* Removes the obsolete `solvers.jl`.


## 0.11.0

### Breaking Changes

* Rename `AtomicSolution` to `SolutionStep`
* Disable `Simulation` functionality until `EnsembleSolution` is added to [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl)

### New Features


### Documentation

* Include documentation of [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl) and [GeometricSolutions.jl](https://github.com/JuliaGNI/GeometricSolutions.jl)


## 0.10.0

### Breaking Changes

* Refactor `TimeSeries`, `DataSeries` and `Solution` and move to [GeometricSolutions.jl](https://github.com/JuliaGNI/GeometricSolutions.jl)
* Adapt Solution HDF5 interface to default Julia argument order and naming conventions
* Extract HDF5 functionality from Solutions into separate data structure
* Remove parallel Solution types


## 0.9.0

### Breaking Changes

* Move `HermiteInterpolation` to Integrators and remove `Interpolation` sub-package
* Move `Equations` submodule to GeometricEquations.jl
* Move `Common`, `Config` and `Utils` submodules to [GeometricBase.jl](https://github.com/JuliaGNI/GeometricBase.jl)
* Move `TimeSeries`, `DataSeries` and `Solution` from `Solutions` types to [GeometricBase.jl](https://github.com/JuliaGNI/GeometricBase.jl)
* Remove parallel DataSeries and Solution types

### New Features

* Implement first and second order Degenerate Variational Integrators (DVIs)
* Add tests for extrapolation methods

### Fixes

* Bugfixes in implicit equations
* Bugfixes in extrapolation methods
* Bugfixes in initial guesses
* Bugfixes in VPRK and VSPARK initialisation
* Bugfixes in `TimeSeries` `getindex` methods

### Documentation

* Add missing docstrings in various places and remove superficial docstrings


## 0.8.0

### Breaking Changes

* Use RungeKutta.jl for most tableaus and coefficients
* Move stochastic integrators to separate package
* Rewrite of most equation types
* Rename `VODE` and `VDAE` to `LODE` and `LDAE` for consistency with `HODE` and `HDAE`
* Add optional fields for the secondary constraint to all *DAE equations

### New Features

* Allow for arbitrary data structures as states (still experimental and not fully supported)
* Add `convert` methods for `PODE` and `HODE` to `ODE` and `SODE`

### Fixes

* Countless minor bugfixes

### Documentation

* Add theoretical background for variational integrators, Runge-Kutta and splitting methods
* Add references for most methods


## 0.7.0

* Use CompactBasisFunctions.jl instead of BasisFunctions submodule
* Use QuadratureRules.jl instead of Quadratures submodule
* Use SimpleSolvers.jl instead of Solvers submodule
* Use GeometricProblems.jl instead of TestProblems submodule


## 0.6.2

* Bugfix release


## 0.6.1

* Bugfix release


## 0.6.0

### Breaking Changes

* Revise tableaus: align constructor names with RungeKutta.jl

### New Features

* Add new Runge-Kutta tableaus
* Generalise Lobatto and Radau tableaus to arbitrary number of stages
* Extend documentation on integrators and tableaus


## 0.5.1

* Update documentation
* Fix HDF5 v0.14 deprecations


## 0.5.0

* Moved repository to JuliaGNI
* Moved CI from Travis to GitHub

### Breaking Changes

* Functions for initial guesses are now called v̄ and f̄ and can be prescribed separately from v and f in PDAE, HDAE, etc.
* Rename SPARK tableau constructors and unify distinct constructors for Lobatto tableaus with different number of stages

### New Features

* Implement SPARK integrator for index-two DAEs
* Implement infrastructure for storing internal variables and solver output to atomic solutions
* Store internal variables of SPARK and VPRK integrators in atomic solution
* Add various five-stage Lobatto tableaus
* Add and clean up SPARK tableaus and add docstrings
* Add functions for checking symplecticity conditions of SPARK tableaus
* Add Aqua.jl tests

### Fixes

* Fix initial guess warnings in tests by prescribing proper functions for v̄ and f̄ in example problems
* Fix update_multiplier() method for SPARK integrators


## 0.4.1

### New Features

* Atomic solutions can now store a NamedTuple of internal variables of the integrator, including nonlinear solver output
* Output of internal variables has been added to VPRK integrators
* Add Gauss-Legendre tableaus for implicit partitioned Runge-Kutta methods

### Fixes

* Revision of integrator type hierarchy


## 0.4.0

### New Integrators

* Runge-Kutta integrators for implicit ODEs (`FIRKimplicit` and `SRKimplicit`)
* Variational Partitioned Runge-Kutta integrator with projection based on internal stages

### Fixes

* Computation of initial guess in *all* implicit integrators


## Open Issues

The migration checklist in `todo.md` is fully worked off — all thirteen items are closed and the
two trailing sections are post-mortems rather than tasks. What follows is what is genuinely open,
taken from the verification report (`docs/src/audit.md`) and the test suite.

### The solver status is available but not acted on

0.18.2 routes every solve through `solve_with_status!` and hands the status to
`check_solver_status`, whose default returns it and does nothing else. That was the deliberate
choice — SimpleSolvers stays the single reporting voice, so no run changes what it prints — but it
means a step that did not converge is still only *reported*, never *acted on*, and the trajectory
continues past the point where it stopped meaning anything with nothing in `sol` to mark it.

The place to act is GeometricIntegratorsBase's `integrate!`, which already handles two of the three
ways a step can go wrong (a `NonlinearSolverException`, and NaNs in the iterate) by warning with
the time step and returning what was computed so far. This is tracked in that package's
`## Open Issues` rather than here, since the hook and the loop both live there; it is named here
because this package is where the consequences would show — see the `@test_broken` methods below,
several of which reach `max_iterations` on every step.

### SLRK is not symplectic (audit finding S15)

The manuscript's proof leaves an uncontrolled term in the multiplier block: the null-vector
condition that kills the velocity term has no counterpart for Λ, and the proof defers to a label
that is never defined. Measured with PoincareInvariants.jl, the first invariant drifts secularly
at O(h^(p+1)) per step, growing linearly with the step count, while a genuinely symplectic
variational integrator holds it to round-off.

Within this ansatz the defect cannot be removed: exact symplecticity, the primary constraint at
every stage, and the constraint at the solution cannot be had together. The defect is gauge
invariant; solvability is not — the two Lobatto IIIA-IIIB pairs go singular on a gauge-equivalent
one-form with a vanishing component.

The family remains useful — it preserves the constraints — but it should not be described as
symplectic.

### The VSPARK projection methods are not symplectic either (S17)

On the massless charged particle, every one of the thirty Gauss-inner methods (ten projections ×
s = 1,2,3) shows a clean O(h^k) defect with k ∈ {3,4} and 100-step drifts up to 1.2e-05. Nothing
is at round-off. Which tableau condition fails is per-projection, not universal.

A sufficient criterion for a *false* pass is established: if every component of ϑ is at most
linear in q, an unprojected variational Runge-Kutta method already lands on the constraint
manifold, so the projection is inert and any projection built on such an inner method comes out
symplectic whichever condition it violates. That explains `PointVorticesLinear`.

**Still unexplained:** why `LotkaVolterra2d` and the two singular gauges nevertheless come out at
round-off. None of them is in the linear-ϑ class — each has a nonlinear first component — and the
projective multipliers are measured between 2.4e-03 and 2.6e+01, so the projection is
demonstrably active. Only the negative half is established.

What a proof would need is the two-form condition Σᵢ b₄ᵢ dΦ̃ᵢ ∧ dΛ̃ᵢ = 0 rather than the pointwise
Φ̃ᵢ = 0; measurement shows neither factor is the mechanism.

### The midpoint projection is not symplectic

Its proof needs R(∞) = +1, while requiring the midpoint to be an internal stage forces
R(∞) = −1, since a midpoint stage means e_mᵀA = bᵀ/2 and hence bᵀA⁻¹e = 2. With the sign that
makes the method solvable, the defect is O(h³); dropping the R(∞) factor instead makes the step
equations unsolvable, the irreducible residual being the midpoint/trapezoidal discrepancy. For a
tableau without a midpoint stage the equations are solvable but the defect is O(h³) as well, so
that sign does not rescue the method either.

Two consequences for the code: `VPRKpSymplectic` is a post projection whose R(∞) is absorbed by
the multiplier and is therefore **the same map as `VPRKpStandard`**; and the mixed term of the
modified two-form as printed in the source manuscript has the wrong sign (with the sign corrected
the form is preserved to round-off).

Note also that both Lotka-Volterra models are unusable for testing symplecticity here: their
multiplier has a single nonzero component whose ϑ component is affine, so every projection
method, including the standard one, comes out exactly symplectic on them.

### `VPRKpInternal` and `VPRKpSecondary` do not run

Both are exported, and both appear in `test/integrators/test_show.jl` (which never solves), but
neither has an `initial_guess!` method for its state layout. Their tests are commented out at
`test/projections/projections_vprk_tests.jl:109` and `:127`. The verdict recorded above on the
internal projection is therefore one about the method, not about executable code.

### `VSPARK(SPARKLobattoIIIBIIIA(2))` has a singular stage system

Restated in 0.18.2, having been recorded here since 0.17.0 as a case that "stalls". It does stall —
the solve stagnates after 3 iterations at rf_a = 3.57e-5 against the `f_abstol = 5.77e-15` it was
asked for, and where it returns it is the one stagnation warning a full test run still prints — but
the stall is a symptom. The stage system is **numerically singular**: cond ≈ 5.6E16, σmin = 5.7E-17
against σmax = 3.2, with σmin an order below the `n·eps·σmax ≈ 1.8E-14` at which a 26×26 system
stops having a numerical rank. That is the same matrix, to within a factor of 1.5 in σmin, as the
`SPARKLobattoIIIAIIIB(2)` and `SPARKGLRK(2)` siblings the suite asserts `SingularException` for.

The 0.17.0 promotion from `@test_broken` to `@test` is therefore retracted: it read one platform's
luck as a property of the method. Whether LAPACK's `getrf` lands on an exact zero pivot or on one
of ~1E-17 decides the outcome, so the same call raises on Julia 1.13 and nightly under Linux and
Windows and returns a ~1E-6 answer on 1.10 and 1.12 everywhere and on macOS throughout. The test
now accepts both. What stays open is the method: `s = 2` is audit finding S8, degenerate at the
lowest stage count, and the answer it sometimes returns is one produced by a Newton direction
solved out of a rank-deficient matrix. It should not be treated as a supported configuration.

### A rank-deficient stage system is diagnosed by luck rather than by design

The entry above is one case of a general gap, and the `Known-broken cases remain` entry below
already names the class — the *marginally* singular methods "which converge or zero-pivot depending
on rounding, so that no single assertion is reliable for them". What makes them unassertable is
that nothing on the path distinguishes rank deficiency from a hard problem. SimpleSolvers' `LU`
linear solver raises `SingularException` only on an exact zero pivot; a pivot of 1E-17 is accepted,
and the Newton direction that comes back out of it is arbitrary in magnitude and direction. The
integrator sees a solve that stalls, not a matrix that has no rank, and the two call for opposite
responses.

A rank-revealing factorization, or simply a pivot-magnitude threshold relative to `σmax`, would
make the outcome deterministic and let a caller distinguish "this method is degenerate here" from
"this step is hard". Both belong in SimpleSolvers rather than here, and neither is a compat-bump
change. Until then the SPARK suite has to accept two outcomes for at least one case, and this
package cannot tell a caller which of the two it got.

### The four state-building `solve_with_status!` sites allocate a state per call

0.18.2 leaves DIRK's per-stage loop and the three projection integrators on the state-*building*
form of `solve_with_status!`, which constructs a `NonlinearSolverState` on every call — once per
stage per step for DIRK, once per step for each projection. This is not a regression: the
`solve!(x, s, params)` they replaced went through the same `NonlinearSolverState(x, value(cache(s)))`
convenience path, so nothing got slower. It is simply now visible, and it is the objection
SimpleSolvers 0.12.1's own docstring raises against that form — *"a caller stepping through time
should build one `NonlinearSolver` and one `NonlinearSolverState` and reuse both"*.

Closing it means giving `ProjectionIntegrator` a `solverstate` field and `SingleStageSolvers` one
state per stage, both structural changes to types that GeometricIntegratorsBase and this package
own respectively. Out of a compat bump, and worth doing together rather than one at a time.

### `check_solver_status` cannot tell a caller *which* solve it is being asked about

The hook takes `(status, int)`, which is the right signature for the fifteen methods that solve
once per step. It is thinner than it should be for the other four. DIRK calls it once per stage
with the same `int` every time, so an override cannot tell which of the `s` stage solves failed,
and sees `s` calls per step where every other method produces one. The three projections pass a
`ProjectionIntegrator`, so an override written as GeometricIntegratorsBase's documentation
suggests — on `GeometricIntegrator{<:MyMethod}` — never sees the projection solve at all, only the
inner integrator's.

Neither is wrong as far as it goes, and no test here depends on the distinction. But a caller
overriding the hook to reject a non-converged step gets a coarser instrument than the call sites
could support, and widening the signature is a GeometricIntegratorsBase decision.

### RungeKutta's barred Lobatto tableaus carry 1E-77 where the plain ones carry exact zeros

Noticed while measuring the stage Jacobians above. `TableauLobattoIIIA(s).a` and
`TableauLobattoIIIB(s).a` have an exactly zero first row, as they should; the adjoint variants do
not:

```julia
julia> TableauLobattoIIIB̄(3).a[1,:]
3-element Vector{Float64}:
 -2.1590421387736112e-78
  8.636168555094445e-78
 -2.1590421387736112e-78
```

The magnitude is `eps(BigFloat)` at the default 256-bit precision, so these are rounding residue
from a coefficient solve carried out in `BigFloat` and surviving the conversion to `Float64`. They
reach this package through every SPARK tableau built on the barred pairs — `SPARKLobattoIIIBIIIA(s)`
for `s ≥ 3` shows them in `tableau.p.a`.

Numerically inert: 1E-77 against coefficients of order 1 changes no arithmetic here. What it does
change is that the structural zeros are no longer *detectable* — `iszero`, `count(iszero, …)` and
anything asking "is the first stage explicit?" answer wrongly on these tableaus. Nothing in this
package asks today. The fix is upstream in RungeKutta.jl, where rounding the solve back to exact
zeros costs nothing.

### The PGLRK status-hook override in the test suite is session-global

`test/verification/pglrk_convergence_tests.jl` adds a counting method to
`GeometricIntegratorsBase.check_solver_status` for `GeometricIntegrator{<:PGLRK}`. A method is
global to the session and `runtests.jl` drives its files with `@safetestset` — a fresh module in
the same process — so it also counts every PGLRK integration in `methods_tests.jl`,
`test_show.jl` and `spark_tableaus_tests.jl`, which run after it. Harmless today: it returns its
argument unchanged and nothing there reads the counter. It is recorded because a second counting
override, for another method or another hook, would silently collide with this one, and because
there is no way to scope a method to a file.

### Known-broken cases remain

A full run reports 26 broken assertions, concentrated in SPARK: 10 in the SPARK convergence suite
and 6 in the SPARK integrator suite, with 3 variational, 3 DGVI, 2 PRK, 1 RK and 1
PGLRK/VPRKpTableau making up the rest. (The suite contains 24 `@test_broken` statements; the two
tallies differ because several sit inside loops and others are not reached.)

Most are inherent properties of the methods rather than implementation defects: symplectic plus
constraint-at-solution reduces order or diverges for the Lobatto IIIA-IIIB and IIIB-IIIA
SPARK/HPARK families; R(∞) = (-1)^(s+1) ≠ 1 drops GLVPRK and HPARKGLRK from order 2s to 2 at
s = 2; coinciding tableau pairs give a singular stage system at s = 2. The SPARK cases still left
as `@test_broken`, rather than asserting a specific failure mechanism, are the *marginally*
singular ones, which converge or zero-pivot depending on rounding, so that no single assertion is
reliable for them.

### `order(VPRKGauss(s))` and `order(VPSRK3())` report wrong orders

Inherited from RungeKutta metadata rather than computed here.

### `HSPARKsecondary` remains EXPERIMENTAL

The `BoundsError` and `SingularException` fixed in 0.16.7 got the family as far as the solver, but
a deeper singularity remains in its ω secondary-constraint block.

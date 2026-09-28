# Known issues

Defects found and recorded, not fixed. Each entry gives its kind and its evidence.

### K1 · The solver status is available but not acted on

- location: —
- evidence: 0.18.2 routes every solve through `solve_with_status!` and hands the status to
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
- kind: defect
- found: 2026-08-15; one issue with GeometricIntegratorsBase `K5`

### K2 · SLRK is not symplectic (audit finding S15)

- location: —
- evidence: The manuscript's proof leaves an uncontrolled term in the multiplier block: the null-vector
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
- kind: defect
- found: 2026-08-15

### K3 · The VSPARK projection methods are not symplectic either (S17)

- location: —
- evidence: On the massless charged particle, every one of the thirty Gauss-inner methods (ten projections ×
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
- kind: defect
- found: 2026-08-15

### K4 · The midpoint projection is not symplectic

- location: —
- evidence: Its proof needs R(∞) = +1, while requiring the midpoint to be an internal stage forces
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
- kind: defect
- found: 2026-08-15

### K5 · `VPRKpInternal` and `VPRKpSecondary` do not run

- location: `test/integrators/test_show.jl`
- evidence: Both are exported, and both appear in `test/integrators/test_show.jl` (which never solves), but
  neither has an `initial_guess!` method for its state layout. Their tests are commented out at
  `test/projections/projections_vprk_tests.jl:109` and `:127`. The verdict recorded above on the
  internal projection is therefore one about the method, not about executable code.
- kind: defect
- found: 2026-08-15

### K6 · `VSPARK(SPARKLobattoIIIBIIIA(2))` has a singular stage system

- location: —
- evidence: Restated in 0.18.2, having been recorded here since 0.17.0 as a case that "stalls". It does stall —
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
- kind: defect
- found: 2026-08-16

### K7 · A rank-deficient stage system is diagnosed by luck rather than by design

- location: —
- evidence: The entry above is one case of a general gap, and the `Known-broken cases remain` entry below
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
- kind: upstream
- found: 2026-08-16

### K8 · The four state-building `solve_with_status!` sites allocate a state per call

- location: —
- evidence: 0.18.2 leaves DIRK's per-stage loop and the three projection integrators on the state-*building*
  form of `solve_with_status!`, which constructs a `NonlinearSolverState` on every call — once per
  stage per step for DIRK, once per step for each projection. This is not a regression: the
  `solve!(x, s, params)` they replaced went through the same `NonlinearSolverState(x, value(cache(s)))`
  convenience path, so nothing got slower. It is simply now visible, and it is the objection
  SimpleSolvers 0.12.1's own docstring raises against that form — *"a caller stepping through time
  should build one `NonlinearSolver` and one `NonlinearSolverState` and reuse both"*.

  Closing it means giving `ProjectionIntegrator` a `solverstate` field and `SingleStageSolvers` one
  state per stage, both structural changes to types that GeometricIntegratorsBase and this package
  own respectively. Out of a compat bump, and worth doing together rather than one at a time.
- kind: defect
- found: 2026-08-16

### K9 · `check_solver_status` cannot tell a caller *which* solve it is being asked about

- location: —
- evidence: The hook takes `(status, int)`, which is the right signature for the fifteen methods that solve
  once per step. It is thinner than it should be for the other four. DIRK calls it once per stage
  with the same `int` every time, so an override cannot tell which of the `s` stage solves failed,
  and sees `s` calls per step where every other method produces one. The three projections pass a
  `ProjectionIntegrator`, so an override written as GeometricIntegratorsBase's documentation
  suggests — on `GeometricIntegrator{<:MyMethod}` — never sees the projection solve at all, only the
  inner integrator's.

  Neither is wrong as far as it goes, and no test here depends on the distinction. But a caller
  overriding the hook to reject a non-converged step gets a coarser instrument than the call sites
  could support, and widening the signature is a GeometricIntegratorsBase decision.
- kind: upstream
- found: 2026-08-16

### K10 · RungeKutta's barred Lobatto tableaus carry 1E-77 where the plain ones carry exact zeros

- location: —
- evidence: Noticed while measuring the stage Jacobians above. `TableauLobattoIIIA(s).a` and
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
- kind: upstream
- found: 2026-08-16

### K11 · The PGLRK status-hook override in the test suite is session-global

- location: `test/verification/pglrk_convergence_tests.jl`
- evidence: `test/verification/pglrk_convergence_tests.jl` adds a counting method to
  `GeometricIntegratorsBase.check_solver_status` for `GeometricIntegrator{<:PGLRK}`. A method is
  global to the session and `runtests.jl` drives its files with `@safetestset` — a fresh module in
  the same process — so it also counts every PGLRK integration in `methods_tests.jl`,
  `test_show.jl` and `spark_tableaus_tests.jl`, which run after it. Harmless today: it returns its
  argument unchanged and nothing there reads the counter. It is recorded because a second counting
  override, for another method or another hook, would silently collide with this one, and because
  there is no way to scope a method to a file.
- kind: defect
- found: 2026-08-16

### K12 · Known-broken cases remain

- location: —
- evidence: A full run reports 26 broken assertions, concentrated in SPARK: 10 in the SPARK convergence suite
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
- kind: defect
- found: 2026-08-15

### K13 · `order(VPRKGauss(s))` and `order(VPSRK3())` report wrong orders

- location: —
- evidence: Inherited from RungeKutta metadata rather than computed here.
- kind: upstream
- found: 2026-08-15

### K14 · `HSPARKsecondary` remains EXPERIMENTAL

- location: —
- evidence: The `BoundsError` and `SingularException` fixed in 0.16.7 got the family as far as the solver, but
  a deeper singularity remains in its ω secondary-constraint block.
- kind: defect
- found: 2026-08-15

### K15 · `SPARK` imports several names from a module that does not own them

- location: `src/SPARK.jl:31-34`
- evidence: `check_all_explicit_imports_via_owners(GeometricIntegrators.SPARK)` (ExplicitImports.jl)
  reports six names imported from `GeometricIntegratorsBase` whose owner, by `Base.which`, is a
  different package: `equation`, `equations` and `timestep` are owned by `GeometricBase`,
  `initialize!` and `method` by `SimpleSolvers`, and `problem` by `GeometricEquations`. Run with
  `JULIA_LOAD_PATH="@:@v1.13:@stdlib" julia --project=<checkout>` and `using ExplicitImports` at
  top level, since ExplicitImports lives in the shared `@v1.13` environment.
- kind: defect
- found: 2026-08-31

### K16 · `src/spark/integrators_spark_parameters.jl` is not included

- location: `src/spark/integrators_spark_parameters.jl`
- evidence: `SPARK.jl` does not `include` this file (`grep -n "include(" src/SPARK.jl`), and no
  other file in `src/` or `test/` names it (`grep -rn integrators_spark_parameters src test`). It
  holds the only definitions of `equation(int::AbstractIntegratorSPARK, i::Symbol)` and
  `equations(int::AbstractIntegratorSPARK)`, both of which `SPARK.jl` imports.
- kind: dead code
- found: 2026-09-27

## KI-1 · `test/helpers/test_functions.jl` is included by no file

- **Kind:** dead code.
- **Evidence:** `grep -rn test_functions test` finds no `include`. On `origin/main`, as
  `test/solutions/test_functions.jl`, no file included it either.

## KI-2 · `splitting_methods_tests.jl` checks the order of two methods only

- **Kind:** missing test.
- **Evidence:** the mutants `order(LieA)` 1→2 and `order(McLachlan2)` 2→3 in
  `src/integrators/splitting/splitting_methods.jl` survive
  `test/integrators/splitting/splitting_methods_tests.jl`, which asserts `order` only for
  `Yoshida6` and `Yoshida8`.

## KI-3 · Aqua's piracy check does not see the `Integrators` submodule

- **Kind:** missing test.
- **Evidence:** the mutant `Base.length(::Symbol) = 0` in `src/integrators/vi/vprk_methods.jl`
  survives `test/quality/aqua.jl`. `Aqua.Piracy.hunt(GeometricIntegrators)` lists the same method
  when it is defined at the top level of `GeometricIntegrators`, and not when it is defined in
  `GeometricIntegrators.Integrators`.

## KI-4 · `test/simulations/simulations_tests.jl` holds no active test

- **Kind:** dead test.
- **Evidence:** every `@test` in the file is commented out; `run-tests.jl` reports 0 tests.

## KI-5 · A stale path of a moved test file under `src/`

- **Kind:** docs.
- **Evidence:** `src/integrators/rk/integrators_pglrk.jl:31` names
  `test/methods/pglrk_coefficients_tests.jl`, which is now
  `test/integrators/rk/pglrk_coefficients_tests.jl`. The test migration changes nothing under
  `src/`.

## KI-6 · Integrator tests do not mirror the subdirectories of `src/integrators/`

- **Kind:** layout.
- **Evidence:** `rk_integrators_tests.jl`, `rk_implicit_integrators_tests.jl`,
  `splitting_integrators_tests.jl`, `variational_integrators_tests.jl`,
  `hamilton_pontryagin_integrators_tests.jl`, `dvi_integrators_tests.jl` and
  `galerkin_integrators_tests.jl` stay in `test/integrators/`, while their sources are in
  `src/integrators/rk/`, `splitting/`, `vi/`, `hpi/`, `dvi/` and `cgvi/`. `test-layout.jl --check`
  checks only that the directory exists in `src/`, so it does not report this.

### K17 · The docstring of `linear_solver_defaults` is on no docs page

- location: `src/integrators/solver_defaults.jl:1`
- evidence: no `Pages` list in `docs/src/modules/integrators.md` names
  `integrators/solver_defaults.jl` (`grep -rn solver_defaults docs/src` finds nothing). The
  default checkdocs reports the docstring as missing; `docs/make.jl:26` keeps `:missing_docs`
  warn-only, so the build warns and does not fail.
- kind: docs
- found: 2026-09-28, critic round 1 of the `LU()` default

### K18 · The `LU()` default is not measured against `LapackLU` at this package's sizes

- location: `src/integrators/solver_defaults.jl:14`
- evidence: the SimpleSolvers 0.14.0 docstring at `linear/linear_solvers.jl:260–265` gives `LU()`
  2× slower than `LapackLU` at n = 64 and 32× slower at n = 768. No measurement of the stage-system
  sizes of this package exists.
- kind: not verified
- found: 2026-09-28, critic round 1 of the `LU()` default

### K19 · The root walk of `solver_defaults.jl` does not check that the default is applied

- location: `test/integrators/solver_defaults.jl:63`
- evidence: the walk checks only that `which(initsolver, …)` for each method root is defined in
  GeometricIntegrators. An override there that leaves out `linear_solver_defaults` passes it; only
  the DIRK and StandardProjection overrides have read-backs of their own.
- kind: missing test
- found: 2026-09-28, critic round 1 of the `LU()` default

### K20 · `integrate` of a `BigFloat` problem fails in GeometricIntegratorsBase

- location: GeometricIntegratorsBase `src/integrate.jl:113`
- evidence: an ODE with a `BigFloat` timespan gives `TypeError: in Type, in parameter, expected
  Int64, got a value of type BigInt` (critic probe `round-1b/probe_defaults.jl`, earlier version).
- kind: upstream
- found: 2026-09-28, critic round 1 of the `LU()` default

### K21 · `test/integrators/solver_defaults.jl` is not run on the Julia 1.10 floor locally

- location: `test/integrators/solver_defaults.jl`
- evidence: `julia +1.10 … run-tests.jl <worktree> integrators/solver_defaults.jl` stops with
  "Could not locate the source code for the StyledStrings package", from a gitignored 1.13
  `Manifest.toml` in the worktree. The `min` CI job is the check.
- kind: not verified
- found: 2026-09-28, critic round 1 of the `LU()` default

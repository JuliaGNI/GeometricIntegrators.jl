# TODO

# Migration: drop the SimpleSolvers warning workarounds

SimpleSolvers 0.9.3/0.10.0 (JuliaGNI/SimpleSolvers.jl#171) removes the cause of the solver and
line-search warning floods, so the log-suppression scaffolding in this repo can largely go. Two
changes matter here:

1. **A line search now shares its solver's `Options`.** Previously every solver built its
   `Linesearch` with an `Options` constructed from nothing but defaults, so `verbosity = 0` never
   reached the line search. `config(linesearch(s)) === config(s)` now holds. This is exactly the
   blocker recorded at `test/spark/spark_integrators_tests.jl:34` — *"keeps its own Options and
   never sees a `verbosity` kwarg passed through integrate"*.
2. **The benign outcomes no longer warn at the default verbosity.** Reaching the merit's round-off
   floor is the *expected* final state of a converged solve, so it is reported only at
   `verbosity ≥ 2`. Measured on GeometricProblems, 12 of 15 integrations went from warning to
   silent with no downstream change at all.

Note the distinction that runs through the whole list: a *converging* solve should now be silent,
while a *genuinely divergent* one still warns — correctly, and once rather than thousands of times.
Suppression is only still appropriate for the second kind.

## Remove the `muffle` scaffolding — 42 call sites, 6 duplicate definitions

- [x] `test/spark/spark_integrators_tests.jl` — 21 sites, `muffle` defined at line 35
      (`NullLogger`). These wrap `@test_broken` cases that are *genuinely* divergent, so they will
      still warn. Prefer passing `verbosity = 0` through `integrate`'s `options` over a logger:
      it silences only SimpleSolvers, whereas `NullLogger` also swallows unrelated failures. Then
      delete `muffle`. Also drop the now-stale comment at line 34.
- [x] `test/projections/projections_vprk_tests.jl` — 6 sites, definition at line 189
      (`SimpleLogger(devnull, Error)`). The comment there already notes *"`NullLogger` would also
      swallow genuine solver failures"* — that concern disappears with `verbosity = 0`.
- [x] `test/verification/pglrk_convergence_tests.jl` — 3 sites, definition at line 90.
- [x] `test/verification/spark_convergence_tests.jl` — 2 sites, definition at line 36.
- [x] `scripts/slrk_verification.jl` — 6 sites, definition at line 64.
- [x] `scripts/vspark_projection_symplecticity.jl` — 4 sites, definition at line 193.

If any suppression survives, **extract one shared helper** instead of the current six copies in
two variants (`NullLogger` in four files, `SimpleLogger(devnull, Error)` in two). The split was
not deliberate.

## Reconsider the tolerances, now that the solver reports the achievable floor

An `f_abstol` below the residual's own round-off floor is still unsatisfiable — that is arithmetic,
not a bug, and 0.9.3 does not change it. What changed is that the solver now *stops* after
`max_stalls = 2` stalled steps and reports the residual it actually achieved against the tolerance
that was requested, instead of spinning to `max_iterations`. Choosing a value is no longer guesswork.

- [ ] `src/spark/abstract.jl:18-20` — `x_suctol=2eps(), f_abstol=8eps(), f_suctol=2eps()`. Run one
      SPARK case and read the achieved `rfₐ` out of the stagnation message, then set `f_abstol`
      an order of magnitude above it.
- [ ] `test/verification/variational_convergence_tests.jl:31,33` and
      `test/verification/hpi_convergence_tests.jl:38,41` — the `f_abstol = 4e-15` relaxations. Check
      whether they are still needed; if so, the comments should cite the measured floor.
- [x] `src/integrators/rk/integrators_pglrk.jl:364` already passes `verbosity=0` — that one is fine
      and can stay as the model for the rest.

## Revise the audit write-up

- [x] `docs/src/audit.md` — the "Third pass" section records the flood and, at line 496, an
      *"Update (SimpleSolvers 0.9.2)"* paragraph. Both are superseded: the API blocker is gone, and
      the remaining warnings are one actionable stagnation message per stuck solve rather than
      thousands. The measured counts in that section should be re-taken or marked historical.

## Now-viable assertions

- [x] `@test_nowarn` was abandoned across these suites because *which* orbits tripped a warning
      depended on platform floating-point details. With `verbosity = 0` reaching the line search
      that is no longer true, so `@test_nowarn` (or a `TestLogger` filtered on
      `nameof(_module) === :SimpleSolvers`) becomes a usable regression tripwire again.
- [ ] `status(solver, state)`, `isconverged` and `isstalled` are now exported, so a `@test_broken`
      case can assert *why* it fails (stagnation at the residual floor) instead of measuring an
      error norm. `solverstate(int)` is already retained, so nothing needs threading through.

## Cross-repo dependency

- [ ] `GeometricIntegratorsBase/src/solvers.jl:6` sets the framework default `f_abstol = 8eps()`,
      which is the root of the original flood and is below the round-off floor of most residuals
      here. Line 8 has `# verbosity=2,` commented out in the same `default_options`. Changing that
      default is a decision for that repo, but it is what would make the *library-default* paths
      quiet — and note that `GeometricIntegratorsBase` replaces its whole `default_options` set as
      soon as a caller passes any option, so callers must restate `min_iterations` etc.

# Known issues

Defects found and recorded, not fixed. Each entry gives its kind and its evidence.

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

## KI-5 · Stale paths of moved test files

- **Kind:** docs.
- **Evidence:** `test/integrators/rk_integrators_tests.jl:22` names
  `test/methods/runge_kutta_methods_tests.jl`, and `test/spark/spark_tableaus_tests.jl:93` and
  `src/integrators/rk/integrators_pglrk.jl:31` name `test/methods/pglrk_coefficients_tests.jl`.
  Both files are now under `test/integrators/rk/`.

## KI-6 · Integrator tests do not mirror the subdirectories of `src/integrators/`

- **Kind:** layout.
- **Evidence:** `rk_integrators_tests.jl`, `rk_implicit_integrators_tests.jl`,
  `splitting_integrators_tests.jl`, `variational_integrators_tests.jl`,
  `hamilton_pontryagin_integrators_tests.jl`, `dvi_integrators_tests.jl` and
  `galerkin_integrators_tests.jl` stay in `test/integrators/`, while their sources are in
  `src/integrators/rk/`, `splitting/`, `vi/`, `hpi/`, `dvi/` and `cgvi/`. `test-layout.jl --check`
  checks only that the directory exists in `src/`, so it does not report this.

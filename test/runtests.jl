using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Method list" include("integrators/method_list_tests.jl")
    @safetestset "Solver defaults" include("integrators/solver_defaults.jl")
    @safetestset "Runge-Kutta methods" include("integrators/rk/runge_kutta_methods_tests.jl")
    @safetestset "PGLRK coefficients" include("integrators/rk/pglrk_coefficients_tests.jl")
    @safetestset "Splitting methods" include("integrators/splitting/splitting_methods_tests.jl")
    @safetestset "VPRK methods" include("integrators/vi/vprk_methods_tests.jl")
    @safetestset "Runge-Kutta integrators for implicit equations" include("integrators/rk_implicit_integrators_tests.jl")
    @safetestset "Splitting integrators" include("integrators/splitting_integrators_tests.jl")
    @safetestset "Degenerate variational integrators" include("integrators/dvi_integrators_tests.jl")
    @safetestset "Hamilton-Pontryagin integrators" include("integrators/hamilton_pontryagin_integrators_tests.jl")
    @safetestset "Show methods and integrators" include("integrators/test_show.jl")
    @safetestset "Ensemble integrators" include("integrators/ensemble_integrators_tests.jl")
    @safetestset "Projection methods" include("projections/projections_tests.jl")
    @safetestset "Projection methods with implicit equations" include("projections/projections_implicit_tests.jl")
    @safetestset "Simulations" include("simulations/simulations_tests.jl")
    @safetestset "SPARK tableaus" include("spark/spark_tableaus_tests.jl")
    @safetestset "Convergence: Runge-Kutta" include("verification/rk_convergence_tests.jl")
    @safetestset "Convergence: partitioned Runge-Kutta" include("verification/prk_convergence_tests.jl")
    @safetestset "Convergence: splitting and composition" include("verification/splitting_convergence_tests.jl")
    @safetestset "Convergence: formal Lagrangian Runge-Kutta" include("verification/flrk_convergence_tests.jl")
    @safetestset "Convergence: projected Gauss-Legendre RK and VPRKpTableau" include("verification/pglrk_convergence_tests.jl")
    @safetestset "Convergence: Galerkin variational integrators" include("verification/galerkin_convergence_tests.jl")
    @safetestset "Convergence: Hamilton-Pontryagin integrators" include("verification/hpi_convergence_tests.jl")
end
if "slow" in GROUPS
    @safetestset "Runge-Kutta integrators" include("integrators/rk_integrators_tests.jl")
    @safetestset "Variational integrators" include("integrators/variational_integrators_tests.jl")
    @safetestset "Galerkin variational integrators" include("integrators/galerkin_integrators_tests.jl")
    @safetestset "Projection methods with VPRK integrators" include("projections/projections_vprk_tests.jl")
    @safetestset "SPARK integrators" include("spark/spark_integrators_tests.jl")
    @safetestset "Convergence: variational integrators" include("verification/variational_convergence_tests.jl")
    @safetestset "Convergence: degenerate variational integrators" include("verification/dvi_convergence_tests.jl")
    @safetestset "Convergence: SPARK integrators" include("verification/spark_convergence_tests.jl")
    @safetestset "Convergence: discontinuous Galerkin variational integrators" include("verification/dgvi_convergence_tests.jl")
    @safetestset "Convergence: projected integrators" include("verification/projection_convergence_tests.jl")
end
if "broken" in GROUPS
    @safetestset "Solution steps" include("integration/solution_step_tests.jl")  # issue #248
    @safetestset "Deterministic solutions" include("integration/deterministic_solutions_tests.jl")  # issue #252
    @safetestset "Solution I/O" include("integration/io_tests.jl")  # issue #253
    @safetestset "Initial guesses (harmonic oscillator)" include("integration/initial_guess_tests_harmonic_oscillator.jl")  # issue #249
    @safetestset "Initial guesses (Lotka-Volterra)" include("integration/initial_guess_tests_lotka_volterra.jl")  # issue #250
end

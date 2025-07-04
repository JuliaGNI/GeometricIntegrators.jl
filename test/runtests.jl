using SafeTestsets

@safetestset "Extrapolation Methods                                                           " begin
    include("extrapolation/extrapolation_tests.jl")
end
@safetestset "Solution Tests                                                                  " begin
    include("solutions/solutions_tests.jl")
end
@safetestset "Initial Guesses (Harmonic Oscillator)                                           " begin
    include("extrapolation/initial_guess_tests_harmonic_oscillator.jl")
end
@safetestset "Initial Guesses (Lotka-Volterra)                                                " begin
    include("extrapolation/initial_guess_tests_lotka_volterra.jl")
end

@safetestset "Euler Integrators                                                               " begin
    include("integrators/euler_tests.jl")
end
@safetestset "Runge-Kutta Integrators                                                         " begin
    include("integrators/rk_integrators_tests.jl")
end
@safetestset "Runge-Kutta Integrators for Implicit Equations                                  " begin
    include("integrators/rk_implicit_integrators_tests.jl")
end
@safetestset "Splitting Integrators                                                           " begin
    include("integrators/splitting_integrators_tests.jl")
end

@safetestset "Variational Integrators                                                         " begin
    include("integrators/variational_integrators_tests.jl")
end
@safetestset "Degenerate Variational Integrators                                              " begin
    include("integrators/dvi_integrators_tests.jl")
end
@safetestset "Galerkin Variational Integrators                                                " begin
    include("integrators/galerkin_integrators_tests.jl")
end
@safetestset "Hamilton-Pontryagin Integrators                                                 " begin
    include("integrators/hamilton_pontryagin_integrators_tests.jl")
end

@safetestset "Projection Methods                                                              " begin
    include("projections/projections_tests.jl")
end
@safetestset "Projection Methods with Implicit Equations                                      " begin
    include("projections/projections_implicit_tests.jl")
end
@safetestset "Projection Methods with Variational Partitioned Runge-Kutta Integrators         " begin
    include("projections/projections_vprk_tests.jl")
end

@safetestset "Ensemble Integrator Tests                                                       " begin
    include("integrators/ensemble_integrators_tests.jl")
end

@safetestset "Method Tests                                                                    " begin
    include("methods/methods_tests.jl")
end
@safetestset "Simulation Tests                                                                " begin
    include("simulations/simulations_tests.jl")
end

@safetestset "SPARK Integrators                                                               " begin
    include("spark/spark_integrators_tests.jl")
end
@safetestset "SPARK Tableau Tests                                                             " begin
    include("spark/spark_tableaus_tests.jl")
end


using SafeTestsets

@safetestset "Extrapolation Methods                                                           " begin include("extrapolation_tests.jl") end
@safetestset "Solution Tests                                                                  " begin include("solutions/solutions_tests.jl") end
@safetestset "Initial Guesses                                                                 " begin include("integrators/initial_guess_tests.jl") end
@safetestset "Runge-Kutta Integrators                                                         " begin include("integrators/rk_integrators_tests.jl") end
@safetestset "Runge-Kutta Integrators for Implicit Equations                                  " begin include("integrators/rk_implicit_integrators_tests.jl") end
@safetestset "Variational Integrators                                                         " begin include("integrators/variational_integrators_tests.jl") end
@safetestset "Degenerate Variational Integrators                                              " begin include("integrators/dvi_integrators_tests.jl") end

# @safetestset "Galerkin Variational Integrators                                                " begin include("integrators/galerkin_integrators_tests.jl") end
# @safetestset "Splitting Integrators                                                           " begin include("integrators/splitting_integrators_tests.jl") end
@safetestset "SPARK Integrators                                                               " begin include("integrators/spark_integrators_tests.jl") end

@safetestset "Projection Methods                                                              " begin include("projections/projections_tests.jl") end
@safetestset "Projection Methods with Implicit Equations                                      " begin include("projections/projections_implicit_tests.jl") end
@safetestset "Projection Methods with Variational Partitioned Runge-Kutta Integrators         " begin include("projections/projections_vprk_tests.jl") end

@safetestset "Common Integrator Functionality                                                 " begin include("integrators/integrators_common_tests.jl") end
@safetestset "Tableau Tests                                                                   " begin include("tableaus/tableaus_tests.jl") end
# @safetestset "Method Tests                                                                    " begin include("methods/methods_tests.jl") end

# @safetestset "Simulation Tests                                                                " begin include("simulations/simulations_tests.jl") end

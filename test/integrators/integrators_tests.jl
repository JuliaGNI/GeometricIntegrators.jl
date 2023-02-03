
using SafeTestsets

# @safetestset "Common Integrator Functionality                                                 " begin include("integrators_common_tests.jl") end
# @safetestset "Initial Guesses                                                                 " begin include("initial_guess_tests.jl") end
# @safetestset "Runge-Kutta Integrators                                                         " begin include("rk_integrators_tests.jl") end
# @safetestset "VPRK Integrators                                                                " begin include("vprk_integrators_tests.jl") end
# @safetestset "SPARK Integrators                                                               " begin include("spark_integrators_tests.jl") end
# @safetestset "Degenerate Variational Integrators                                              " begin include("dvi_integrators_tests.jl") end
# @safetestset "Special Integrators                                                             " begin include("special_integrators_tests.jl") end
# @safetestset "Splitting Integrators                                                           " begin include("splitting_integrators_tests.jl") end
# @safetestset "Galerkin Integrators                                                            " begin include("galerkin_integrators_tests.jl") end

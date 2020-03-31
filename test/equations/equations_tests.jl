
using SafeTestsets

@safetestset "General Equation Functionality                                                  " begin include("general_equations_tests.jl") end
@safetestset "Deterministic Equations                                                         " begin include("deterministic_equations_tests.jl") end
@safetestset "Stochastic Equations                                                            " begin include("stochastic_equations_tests.jl") end
@safetestset "Problem Equations                                                               " begin include("problem_equations_tests.jl") end

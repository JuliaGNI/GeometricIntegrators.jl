
using SafeTestsets

@safetestset "Linear Solvers                                                                  " begin include("linear_solvers_tests.jl") end
@safetestset "Jacobians                                                                       " begin include("jacobian_tests.jl") end
@safetestset "Nonlinear Solvers                                                               " begin include("nonlinear_solvers_tests.jl") end


using SafeTestsets

@safetestset "Solver Tests                                                                    " begin include("solvers/solvers_tests.jl") end
@safetestset "Basis Function Tests                                                            " begin include("basis_functions/basis_functions_tests.jl") end
@safetestset "Equation Tests                                                                  " begin include("equations/equations_tests.jl") end
@safetestset "Solution Tests                                                                  " begin include("solution/solution_tests.jl") end
@safetestset "Tableau Tests                                                                   " begin include("tableaus/tableaus_tests.jl") end
@safetestset "Integrator Tests                                                                " begin include("integrators/integrators_tests.jl") end

#include("utils/hdf5_tests.jl")

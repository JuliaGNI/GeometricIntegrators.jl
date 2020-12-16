
using SafeTestsets

@safetestset "Equation Tests                                                                  " begin include("equations/equations_tests.jl") end
@safetestset "Solution Tests                                                                  " begin include("solutions/solutions_tests.jl") end
@safetestset "Tableau Tests                                                                   " begin include("tableaus/tableaus_tests.jl") end
@safetestset "Integrator Tests                                                                " begin include("integrators/integrators_tests.jl") end
@safetestset "Simulation Tests                                                                " begin include("simulations/simulations_tests.jl") end

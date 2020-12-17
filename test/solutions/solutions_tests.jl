
using SafeTestsets

@safetestset "Data Series                                                                     " begin include("dataseries_tests.jl") end
@safetestset "Time Series                                                                     " begin include("timeseries_tests.jl") end
@safetestset "Atomic Solutions                                                                " begin include("atomic_solutions_tests.jl") end
@safetestset "Abstract Solution                                                               " begin include("abstract_solution_tests.jl") end
@safetestset "Deterministic Solutions                                                         " begin include("deterministic_solutions_tests.jl") end
# @safetestset "Stochastic Solutions                                                            " begin include("stochastic_solutions_tests.jl") end
# TODO # reactivate

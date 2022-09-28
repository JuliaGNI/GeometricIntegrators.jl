
using SafeTestsets

@safetestset "Solution Step                                                                   " begin
    include("solution_step_tests.jl")
end
@safetestset "Deterministic Solutions                                                         " begin
    include("deterministic_solutions_tests.jl")
end

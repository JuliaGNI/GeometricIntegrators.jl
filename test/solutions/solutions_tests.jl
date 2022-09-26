
using SafeTestsets

@safetestset "Atomic Solutions                                                                " begin
    include("atomic_solutions_tests.jl")
end
@safetestset "Deterministic Solutions                                                         " begin
    include("deterministic_solutions_tests.jl")
end

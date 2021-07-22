
using GeometricBase
using GeometricIntegrators.Solutions
using Test


struct TestSolution{dType, tType, N} <: Solution{dType, tType, N} end

test_sol = TestSolution{Float64, Float64, 2}()

@testset "$(rpad("Interface Definition",80))" begin
    @test_throws ErrorException timesteps(test_sol)
    @test_throws ErrorException hdf5(test_sol)
    @test_throws ErrorException ntime(test_sol)
    @test_throws ErrorException nsave(test_sol)
    @test_throws ErrorException offset(test_sol)
    @test_throws ErrorException create_hdf5(test_sol, "test.hdf5")
end


using GeometricBase.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators: CacheType, Parameters
using Test


struct MyParameters{DT,TT} <: Parameters{DT,TT} end
struct MyIntegrator{DT,TT} <: AbstractIntegrator{DT,TT} end

params = MyParameters{Float64,Float64}()
int    = MyIntegrator{Float64,Float64}()

# @test_throws ErrorException IntegratorCache(params)
# @test_throws ErrorException IntegratorCache(int)
# @test_throws ErrorException CacheType(Float32, params)

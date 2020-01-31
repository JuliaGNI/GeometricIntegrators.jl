
using GeometricIntegrators
using GeometricIntegrators.Solutions
using Test


@testset "$(rpad("Timeseries",80))" begin
    ntime = 10
    Δt = 0.1

    ts = TimeSeries{eltype(Δt)}(ntime, Δt, 1)
    @test typeof(ts) <: AbstractArray
    @test firstindex(ts)   == 0
    @test firstindex(ts,1) == firstindex(ts.t,1) - 1
    @test firstindex(ts,2) == 1
    @test lastindex(ts)    == ntime
    @test lastindex(ts,1)  == lastindex(ts.t,1) - 1
    @test lastindex(ts,2)  == 1
    @test axes(ts)   == (0:ntime,)
    @test axes(ts,1) == 0:ntime
    @test axes(ts,2) == 1:1
    @test size(ts.t) == (ntime+1,)

    compute_timeseries!(ts, 0.)
    t = collect(0:Δt:ntime*Δt)
    @test ts.t ≈ t atol=eps()

    ts1 = TimeSeries(ntime, Δt)
    compute_timeseries!(ts1, 0.)
    @test ts1 == ts
end

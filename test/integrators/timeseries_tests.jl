
ntime = 10
Δt = 0.1

ts = Timeseries{eltype(Δt)}(ntime, Δt, 1)
@test typeof(ts) <: AbstractArray
@test length(ts) == ntime+1
@test size(ts) == (ntime+1,)
@test indices(ts, 1) == 0:ntime

compute_timeseries!(ts)
t = collect(0:Δt:ntime*Δt)
@test_approx_eq_eps(ts, t, eps())

ts1 = Timeseries(ntime, Δt)
compute_timeseries!(ts1)
@test ts1 == ts

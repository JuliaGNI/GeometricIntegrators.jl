
ntime = 10
Δt = 0.1

ts = TimeSeries{eltype(Δt)}(ntime, Δt, 1)
@test typeof(ts) <: AbstractArray

compute_timeseries!(ts, 0.)
t = collect(0:Δt:ntime*Δt)
@test ts.t ≈ t atol=eps()

ts1 = TimeSeries(ntime, Δt)
compute_timeseries!(ts1, 0.)
@test ts1 == ts

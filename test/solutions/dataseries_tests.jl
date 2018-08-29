
dt = Float64
nd = 2
nt = 10


ni = 1
ds = SDataSeries(dt, nd, nt, ni)
@test typeof(ds) <: AbstractArray{dt,2}
@test lastindex(ds) == (nd, nt)

for i in 0:nt
    ds[1,i] = i
    ds[2,i] = 2^i
end

@test ds.d[1,1] == ds[1,0]
@test ds.d[1:ds.nd,1] == ds[1:ds.nd,0]

@test ds[1,end] == ds.d[1,end]
@test ds[1,nt] == ds[1,end]
@test ds[1:ds.nd,nt] == ds[1:ds.nd,end]

reset!(ds)
@test ds[1,0] == ds[1,end]


ni = 2
ds = SDataSeries(dt, nd, nt, ni)
@test typeof(ds) <: AbstractArray{dt,3}
@test lastindex(ds) == (nd, nt, ni)

for j in 1:ni
    for i in 1:nt
        ds[1,i,j] = j*i
    end
end

@test ds.d[1,1,1] == ds[1,0,1]
@test ds.d[1:ds.nd,1,1] == ds[1:ds.nd,0,1]

@test ds[1,end,1] == ds.d[1,end,1]
@test ds[1,nt,1] == ds[1,end,1]
@test ds[1:ds.nd,nt,1] == ds[1:ds.nd,end,1]

reset!(ds)
@test ds[1,0,1] == ds[1,end,1]

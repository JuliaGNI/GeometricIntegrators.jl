
dt = Float64
nd = 2
nt = 10


ni = 1
ds = SDataSeries(dt, nd, nt, ni)
@test typeof(ds) <: AbstractArray{dt,2}
@test lastindex(ds) == (nd, nt)
@test lastindex(ds,1) == lastindex(ds.d,1)
@test lastindex(ds,2) == lastindex(ds.d,2) - 1
@test axes(ds) == (1:nd, 0:nt)
@test axes(ds,1) == 1:nd
@test axes(ds,2) == 0:nt
@test size(ds.d) == (nd, nt+1)

for i in 0:nt
    ds[1,i] = i
    ds[2,i] = 2^i
end

@test ds.d[1,1] == ds[1,0]
@test ds.d[1:ds.nd,1] == ds[1:ds.nd,0]

@test ds[1,end] == ds.d[1,end]
@test ds[1,nt] == ds[1,end]
@test ds[1:ds.nd,nt] == ds[1:ds.nd,end]

@test ds.d[:,1] == ds[:,0]

reset!(ds)
@test ds[1,0] == ds[1,end]


ni = 2
ds = SDataSeries(dt, nd, nt, ni)
@test typeof(ds) <: AbstractArray{dt,3}
@test lastindex(ds) == (nd, nt, ni)
@test lastindex(ds,1) == lastindex(ds.d,1)
@test lastindex(ds,2) == lastindex(ds.d,2) - 1
@test lastindex(ds,3) == lastindex(ds.d,3)
@test axes(ds) == (1:nd, 0:nt, 1:ni)
@test axes(ds,1) == 1:nd
@test axes(ds,2) == 0:nt
@test axes(ds,3) == 1:ni
@test size(ds.d) == (nd, nt+1, ni)


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

@test ds.d[:,1,1] == ds[:,0,1]
@test ds.d[:,1,:] == ds[:,0,:]

reset!(ds)
@test ds[1,0,1] == ds[1,end,1]

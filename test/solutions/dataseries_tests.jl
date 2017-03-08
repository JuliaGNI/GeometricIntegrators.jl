
dt = Float64
nd = 2
nt = 10


ni = 1
ds = SDataSeries(dt, nd, nt, ni)
@test typeof(ds) <: AbstractArray{dt,2}
@test length(ds) == nd*(nt+1)
@test size(ds) == (nd, nt+1)
@test indices(ds) == (1:nd, 0:nt)

for i in 1:nt
    ds[1,i] = i
end

@test ds.d[1,1] == ds[1,0]
@test ds.d[1:ds.nd,1] == ds[1:ds.nd,0]


ni = 2
ds = SDataSeries(dt, nd, nt, ni)
@test typeof(ds) <: AbstractArray{dt,3}
@test length(ds) == nd*(nt+1)*ni
@test size(ds) == (nd, nt+1, ni)
@test indices(ds) == (1:nd, 0:nt, 1:ni)

for j in 1:ni
    for i in 1:nt
        ds[1,i,j] = j*i
    end
end

@test ds.d[1,1,1] == ds[1,0,1]
@test ds.d[1:ds.nd,1,1] == ds[1:ds.nd,0,1]

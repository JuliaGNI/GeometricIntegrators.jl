
using GeometricIntegrators
using GeometricIntegrators.Solutions
using Test


dt = Float64
nt = 10

@testset "$(rpad("Dataseries 1d",80))" begin
    nd = 1
    ni = 1

    ds = SDataSeries(dt, nt)
    @test typeof(ds) != typeof(SDataSeries(dt, nd, nt))
    @test typeof(ds) != typeof(SDataSeries(dt, nd, nt, ni))
    @test typeof(ds) <: AbstractArray{dt,1}
    @test firstindex(ds)   == 0
    @test firstindex(ds,1) == firstindex(ds.d,1) - 1
    @test firstindex(ds,2) == 1
    @test lastindex(ds)    == nt
    @test lastindex(ds,1)  == lastindex(ds.d,1) - 1
    @test lastindex(ds,2)  == 1
    @test strides(ds)      == (1,)
    @test stride(ds,1)     == 1
    @test axes(ds)   == (0:nt,)
    @test axes(ds,1) == 0:nt
    @test axes(ds,2) == 1:1
    @test size(ds.d) == (nt+1,)
    @test size(ds.d) == size(ds)
    @test size(ds.d,1) == nt+1
    @test ndims(ds)  == 1
    @test eltype(ds) == dt
    @test parent(ds) == ds.d

    @test eachindex(ds) == 0:nt
    @test eachindex(IndexLinear(), ds) == 0:nt
    @test eachindex(IndexCartesian(), ds) == CartesianIndices((0:nt,))

    for i in 0:nt
        ds[i] = i
    end

    @test ds.d[1] == ds[0]
    @test ds.d[end] == ds[nt]
    @test ds.d[end] == ds[end]

    reset!(ds)
    @test ds[0] == ds[end]

    d = rand(nt+1)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)
end


@testset "$(rpad("Dataseries 2d",80))" begin
    nd = 2
    ni = 1

    ds = SDataSeries(dt, nd, nt, ni)
    @test typeof(ds) == typeof(SDataSeries(dt, nd, nt))
    @test typeof(ds) <: AbstractArray{dt,2}
    @test firstindex(ds)   == firstindex(ds.d)
    @test firstindex(ds,1) == firstindex(ds.d,1)
    @test firstindex(ds,2) == firstindex(ds.d,2) - 1
    @test firstindex(ds,3) == 1
    @test lastindex(ds)    == lastindex(ds.d)
    @test lastindex(ds,1)  == lastindex(ds.d,1)
    @test lastindex(ds,2)  == lastindex(ds.d,2) - 1
    @test lastindex(ds,3)  == 1
    @test strides(ds)      == (1,nd)
    @test stride(ds,1)     == 1
    @test axes(ds)   == (1:nd, 0:nt)
    @test axes(ds,1) == 1:nd
    @test axes(ds,2) == 0:nt
    @test axes(ds,3) == 1:1
    @test size(ds.d) == (nd, nt+1)
    @test size(ds.d) == size(ds)
    @test size(ds.d,1) == nd
    @test size(ds.d,2) == nt+1
    @test ndims(ds)  == 2
    @test eltype(ds) == dt
    @test parent(ds) == ds.d

    for i in 0:nt
        ds[1,i] = i
        ds[2,i] = 2^i
    end

    tx = zeros(nd)
    for i in 0:nt
        get_data!(ds, tx, i)
        @test tx == Array{eltype(tx)}([i, 2^i])
    end

    @test ds.d[1,1] == ds[1,0]
    @test ds.d[1:ds.nd,1] == ds[1:ds.nd,0]

    @test ds[1,end] == ds.d[1,end]
    @test ds[1,end] == ds[1,nt]
    @test ds[1:ds.nd,nt] == ds[1:ds.nd,end]

    @test ds.d[:,1] == ds[:,0]

    reset!(ds)
    @test ds[1,0] == ds[1,end]

    tx = rand(nd)
    ds[:,0] .= tx
    @test ds.d[:,1] == tx

    d = rand(nd,nt+1)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)
end


@testset "$(rpad("Dataseries 3d",80))" begin
    nd = 2
    ni = 2

    ds = SDataSeries(dt, nd, nt, ni)
    @test typeof(ds) <: AbstractArray{dt,3}
    @test firstindex(ds)   == firstindex(ds.d)
    @test firstindex(ds,1) == firstindex(ds.d,1)
    @test firstindex(ds,2) == firstindex(ds.d,2) - 1
    @test firstindex(ds,3) == firstindex(ds.d,3)
    @test firstindex(ds,4) == 1
    @test lastindex(ds)    == lastindex(ds.d)
    @test lastindex(ds,1)  == lastindex(ds.d,1)
    @test lastindex(ds,2)  == lastindex(ds.d,2) - 1
    @test lastindex(ds,3)  == lastindex(ds.d,3)
    @test lastindex(ds,4)  == 1
    @test strides(ds)      == (1,nd,nd*(nt+1))
    @test stride(ds,1)     == 1
    @test stride(ds,2)     == nd
    @test axes(ds)   == (1:nd, 0:nt, 1:ni)
    @test axes(ds,1) == 1:nd
    @test axes(ds,2) == 0:nt
    @test axes(ds,3) == 1:ni
    @test axes(ds,4) == 1:1
    @test size(ds.d) == (nd, nt+1, ni)
    @test size(ds.d) == size(ds)
    @test size(ds.d,1) == nd
    @test size(ds.d,2) == nt+1
    @test size(ds.d,3) == ni
    @test ndims(ds)  == 3
    @test eltype(ds) == dt
    @test parent(ds) == ds.d

    for j in 1:ni
        for i in 1:nt
            ds[1,i,j] = j*i
            ds[2,i,j] = j*i^2
        end
    end

    tx = zeros(nd)
    for j in 1:ni
        for i in 1:nt
            get_data!(ds, tx, i, j)
            @test tx == Array{eltype(tx)}([j*i, j*i^2])
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

    tx = rand(nd)
    ds[:,0,1] .= tx
    @test ds.d[:,1,1] == tx

    d = rand(nd,nt+1,ni)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)
end

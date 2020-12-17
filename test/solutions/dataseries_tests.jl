
using GeometricIntegrators
using GeometricIntegrators.Solutions
using Test


dt = Float64
nt = 10

@testset "$(rpad("Dataseries 1d (scalar data)",80))" begin
    nd = 1
    ni = 1

    ds = SDataSeries(rand(dt), nt)
    @test typeof(ds) == typeof(SDataSeries{dt,1}(nt, ni))
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
    @test ds.d[:] == parent(ds[:])
    @test ds.d == collect(0:nt)

    reset!(ds)
    @test ds[0] == ds[end]

    d = rand(nt+1)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)

    for i in 0:nt
        set_data!(ds, dt(i), i)
    end

   @test ds.d == collect(0:nt)
end


@testset "$(rpad("Dataseries 1d (vector-valued data)",80))" begin
    nd = 2
    ni = 1

    ds = SDataSeries(rand(dt, nd), nt)
    @test typeof(ds) == typeof(SDataSeries{Vector{dt},1}(nt, ni))
    @test typeof(ds) <: AbstractArray{Vector{dt},1}
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
    @test eltype(ds) == Vector{dt}
    @test parent(ds) == ds.d

    @test eachindex(ds) == 0:nt
    @test eachindex(IndexLinear(), ds) == 0:nt
    @test eachindex(IndexCartesian(), ds) == CartesianIndices((0:nt,))

    for i in 0:nt
        ds[i] = [i, i^2]
    end

    @test ds.d[1] == ds[0]
    @test ds.d[end] == ds[nt]
    @test ds.d[end] == ds[end]
    @test ds.d[:] == parent(ds[:])
    @test ds.d == [Vector{dt}([i, i^2]) for i in 0:nt]
    
    reset!(ds)
    @test ds[0] == ds[end]

    d = rand(nt+1)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)

    for i in 0:nt
        set_data!(ds, dt(i), i)
    end

   @test ds.d == collect(0:nt)
end


@testset "$(rpad("Dataseries 2d (scalar data)",80))" begin
    ni = 2

    ds = SDataSeries(rand(dt, ni), nt, ni)
    @test typeof(ds) == typeof(SDataSeries{dt,2}(nt, ni))
    @test typeof(ds) <: AbstractArray{dt,2}
    @test firstindex(ds)   == firstindex(ds.d)
    @test firstindex(ds,1) == firstindex(ds.d,1) - 1
    @test firstindex(ds,2) == 1
    @test lastindex(ds)    == lastindex(ds.d)
    @test lastindex(ds,1)  == lastindex(ds.d,1) - 1
    @test lastindex(ds,2)  == lastindex(ds.d,2)
    @test strides(ds)      == (1,nt+1)
    @test stride(ds,1)     == 1
    @test stride(ds,2)     == nt+1
    @test axes(ds)   == (0:nt, 1:ni)
    @test axes(ds,1) == 0:nt
    @test axes(ds,2) == 1:2
    @test size(ds.d) == (nt+1, ni)
    @test size(ds.d) == size(ds)
    @test size(ds.d,1) == nt+1
    @test size(ds.d,2) == ni
    @test ndims(ds)  == 2
    @test eltype(ds) == dt
    @test parent(ds) == ds.d

    for i in 0:nt
        set_data!(ds, dt(i), i, 1)
        set_data!(ds, dt(i)^2, i, 2)
    end

    @test ds.d[:,1] == collect(0:nt)
    @test ds.d[:,2] == collect(0:nt) .^ 2

    @test ds.d[1,1] == ds[0,1]
    @test ds.d[1,1:ni] == ds[0,1:nsamples(ds)]

    @test ds[end,1] == ds.d[end,1]
    @test ds[end,1] == ds[nt,1]
    @test ds[end,1:nsamples(ds)] == ds[nt,1:ni]

    @test ds.d[1,:] == ds[0,:]
    @test ds.d[:,1] == parent(ds[:,1])
    @test ds.d[:,:] == parent(ds[:,:])

    reset!(ds)
    @test ds[0,1] == ds[end,1]

    # tx = rand(nd)
    # ds[:,0] .= tx
    # @test ds.d[:,1] == tx

    d = rand(nt+1,ni)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)
end


@testset "$(rpad("Dataseries 2d (vector-valued data)",80))" begin
    nd = 2
    ni = 5

    ds = SDataSeries([rand(dt, nd) for i in 1:ni], nt)
    @test typeof(ds) == typeof(SDataSeries{Vector{dt},2}(nt, ni))
    @test typeof(ds) <: AbstractArray{Vector{dt},2}
    @test firstindex(ds)   == firstindex(ds.d)
    @test firstindex(ds,1) == firstindex(ds.d,1) - 1
    @test firstindex(ds,2) == 1
    @test lastindex(ds)    == lastindex(ds.d)
    @test lastindex(ds,1)  == lastindex(ds.d,1) - 1
    @test lastindex(ds,2)  == lastindex(ds.d,2)
    @test strides(ds)      == (1,nt+1)
    @test stride(ds,1)     == 1
    @test stride(ds,2)     == nt+1
    @test axes(ds)   == (0:nt, 1:ni)
    @test axes(ds,1) == 0:nt
    @test axes(ds,2) == 1:ni
    @test size(ds.d) == (nt+1, ni)
    @test size(ds.d) == size(ds)
    @test size(ds.d,1) == nt+1
    @test size(ds.d,2) == ni
    @test ndims(ds)  == 2
    @test eltype(ds) == Vector{dt}
    @test parent(ds) == ds.d

    for k in 1:ni
        for i in 0:nt
            ds[i,k][1] = k*i
            ds[i,k][2] = k*i^2
        end
    end

    tx = zeros(dt, nd)
    for k in 1:ni
        for i in 0:nt
            get_data!(ds, tx, i, k)
            @test tx == Array{eltype(tx)}(k .* [i, i^2])
        end
    end

    @test ds.d[1,1] == ds[0,1]
    @test ds.d[1,1:ni] == ds[0,1:nsamples(ds)]

    @test ds[end,1] == ds.d[end,1]
    @test ds[end,1] == ds[nt,1]
    @test ds[end,1:nsamples(ds)] == ds[nt,1:ni]

    @test ds.d[1,:] == ds[0,:]
    @test ds.d[:,1] == parent(ds[:,1])
    @test ds.d[:,:] == parent(ds[:,:])

    reset!(ds)
    @test ds[0,1] == ds[end,1]

    # tx = rand(nd)
    # ds[:,0] .= tx
    # @test ds.d[:,1] == tx

    d = rand(nt+1,ni)
    ds = SDataSeries(d)
    @test ds.d == d

    @test similar(ds).d == zero(d)
end



# TODO: Add tests for array-valued data



# @testset "$(rpad("Dataseries 3d",80))" begin
#     nd = 2
#     ni = 2

#     ds = SDataSeries(dt, nd, nt, ni)
#     @test typeof(ds) <: AbstractArray{dt,3}
#     @test firstindex(ds)   == firstindex(ds.d)
#     @test firstindex(ds,1) == firstindex(ds.d,1)
#     @test firstindex(ds,2) == firstindex(ds.d,2) - 1
#     @test firstindex(ds,3) == firstindex(ds.d,3)
#     @test firstindex(ds,4) == 1
#     @test lastindex(ds)    == lastindex(ds.d)
#     @test lastindex(ds,1)  == lastindex(ds.d,1)
#     @test lastindex(ds,2)  == lastindex(ds.d,2) - 1
#     @test lastindex(ds,3)  == lastindex(ds.d,3)
#     @test lastindex(ds,4)  == 1
#     @test strides(ds)      == (1,nd,nd*(nt+1))
#     @test stride(ds,1)     == 1
#     @test stride(ds,2)     == nd
#     @test axes(ds)   == (1:nd, 0:nt, 1:ni)
#     @test axes(ds,1) == 1:nd
#     @test axes(ds,2) == 0:nt
#     @test axes(ds,3) == 1:ni
#     @test axes(ds,4) == 1:1
#     @test size(ds.d) == (nd, nt+1, ni)
#     @test size(ds.d) == size(ds)
#     @test size(ds.d,1) == nd
#     @test size(ds.d,2) == nt+1
#     @test size(ds.d,3) == ni
#     @test ndims(ds)  == 3
#     @test eltype(ds) == dt
#     @test parent(ds) == ds.d

#     for j in 1:ni
#         for i in 1:nt
#             ds[1,i,j] = j*i
#             ds[2,i,j] = j*i^2
#         end
#     end

#     tx = zeros(nd)
#     ty = zeros(nd,ni)
#     tz = zeros(nd,ni)
#     for i in 1:nt
#         for j in 1:ni
#             get_data!(ds, tx, i, j)
#             @test tx == Array{eltype(tx)}([j*i, j*i^2])
#             tz[1,j] = j*i
#             tz[2,j] = j*i^2
#         end
#         get_data!(ds, ty, i)
#         @test ty == tz
#     end

#     @test ds.d[1,1,1] == ds[1,0,1]
#     @test ds.d[1:ds.nd,1,1] == ds[1:ds.nd,0,1]

#     @test ds[1,end,1] == ds.d[1,end,1]
#     @test ds[1,nt,1] == ds[1,end,1]
#     @test ds[1:ds.nd,nt,1] == ds[1:ds.nd,end,1]

#     @test ds.d[:,1,1] == ds[:,0,1]
#     @test ds.d[:,1,:] == ds[:,0,:]

#     reset!(ds)
#     @test ds[1,0,1] == ds[1,end,1]

#     tx = rand(nd)
#     ds[:,0,1] .= tx
#     @test ds.d[:,1,1] == tx

#     d = rand(nd,nt+1,ni)
#     ds = SDataSeries(d)
#     @test ds.d == d

#     @test similar(ds).d == zero(d)
# end

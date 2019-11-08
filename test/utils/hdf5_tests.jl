
@testset "$(rpad("HDF5 Utils",80))" begin

    using HDF5

    ntime = 10
    nsave = 1
    dim   = 1

    sol = SolutionODE{Float64}(dim, ntime, nsave)

    for i in 1:ntime
        sol[1,i] = i
    end

    tmp = mktempdir()
    file = tmp * "/test.h5"
    writeSolutionToHDF5(sol, file)
    h5 = h5open(file, "r")
    @test h5["x"][:,:] == sol.x
    close(h5)
    rm(tmp, recursive=true)

    tmp = mktempdir()
    h5 = createHDF5(sol, tmp * "/test.h5", 2ntime+1)
    writeSolutionToHDF5(sol, h5)
    writeSolutionToHDF5(sol, h5, ntime)
    close(h5)
    rm(tmp, recursive=true)

end

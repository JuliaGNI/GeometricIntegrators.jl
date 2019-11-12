
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricIntegrators.TestProblems.KuboOscillator
using HDF5: HDF5File
using Test


ntime  = 10
Δt     = .1
h5file = "test.hdf5"


@testset "$(rpad("SDE Solution",80))" begin
    sde  = kubo_oscillator_sde_1()
    ssol = Solution(sde, Δt, ntime)
    @test typeof(ssol) <: SolutionSDE

    # test hdf5 in- and output
    h5 = createHDF5(ssol, h5file)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)

    h5 = createHDF5(ssol, h5file, overwrite=false)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)
    rm(h5file)
end


@testset "$(rpad("PSDE Solution",80))" begin
    psde  = kubo_oscillator_psde_1()
    ssol = Solution(psde, Δt, ntime)
    @test typeof(ssol) <: SolutionPSDE

    # test hdf5 in- and output
    h5 = createHDF5(ssol, h5file)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)

    h5 = createHDF5(ssol, h5file, overwrite=false)
    @test typeof(h5) == HDF5File
    close(h5)
    @test isfile(h5file)
    rm(h5file)
end


@testset "$(rpad("SPSDE Solution",80))" begin
    spsde  = kubo_oscillator_spsde_1()
    ssol = Solution(spsde, Δt, ntime)
    @test typeof(ssol) <: SolutionPSDE
end

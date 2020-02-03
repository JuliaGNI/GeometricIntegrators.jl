
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solutions: createHDF5
using GeometricIntegrators.TestProblems.KuboOscillatorProblem
using HDF5: HDF5File
using Test


nd = 3
ns = 5
nt = 10
Δt = .1

dType  = Float64
h5file = "test.hdf5"


@testset "$(rpad("Wiener Process",80))" begin

    wp = WienerProcess(dType, nd, nt, 1,  Δt, :strong)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :strong)

    wp = WienerProcess(dType, nd, nt, 1,  Δt, :weak)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :weak)

    wp = WienerProcess(dType, nd, nt, ns, Δt, :strong)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :strong)

    wp = WienerProcess(dType, nd, nt, ns, Δt, :weak)
    @test wp == WienerProcess(Δt, wp.ΔW, wp.ΔZ, :weak)

end


@testset "$(rpad("SDE Solution",80))" begin
    sde  = kubo_oscillator_sde_1()
    ssol = Solution(sde, Δt, nt)
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
    ssol = Solution(psde, Δt, nt)
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
    ssol = Solution(spsde, Δt, nt)
    @test typeof(ssol) <: SolutionPSDE
end


using GeometricIntegrators
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using Test

set_config(:nls_atol_break, Inf)
set_config(:nls_rtol_break, Inf)
set_config(:nls_stol_break, Inf)


Δt = 0.1
nt = 100000

nsave1 = 1
nwrte1 = 1000

nsave2 = 10
nwrte2 = 1000

h5file = "test.hdf5"


tab = getTableauImplicitMidpoint()
ode = harmonic_oscillator_ode()


@testset "$(rpad("Serial Simulation",80))" begin

    sim1 = Simulation(ode, tab, Δt, "Harmonic Oscillator Test 1", h5file, nt; nsave=nsave1, nwrite=nwrte1)
    run!(sim1)
    @test isfile(h5file)
    sol1 = SSolutionODE(h5file)
    @test sol1.ntime == nt
    @test sol1.nsave == nsave1
    rm(h5file)

    sim2 = Simulation(ode, tab, Δt, "Harmonic Oscillator Test 2", h5file, nt; nsave=nsave2, nwrite=nwrte2)
    run!(sim2)
    @test isfile(h5file)
    sol2 = SSolutionODE(h5file)
    @test sol2.ntime == nt
    @test sol2.nsave == nsave2
    rm(h5file)

end



nt = 10000
ns = 10

ode = harmonic_oscillator_ode(vcat(rand(1,ns), zeros(1,ns)))


@testset "$(rpad("Parallel Simulation",80))" begin

    sim1 = ParallelSimulation(ode, tab, Δt, "Harmonic Oscillator Test 1", h5file, nt; nsave=nsave1, nwrite=nwrte1)
    run!(sim1)
    @test isfile(h5file)
    sol1 = SSolutionODE(h5file)
    @test sol1.ntime == nt
    @test sol1.nsave == nsave1
    rm(h5file)

    sim2 = ParallelSimulation(ode, tab, Δt, "Harmonic Oscillator Test 2", h5file, nt; nsave=nsave2, nwrite=nwrte2)
    run!(sim2)
    @test isfile(h5file)
    sol2 = SSolutionODE(h5file)
    @test sol2.ntime == nt
    @test sol2.nsave == nsave2
    rm(h5file)

end

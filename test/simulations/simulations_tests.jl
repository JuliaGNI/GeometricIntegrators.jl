
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


Δt = 0.1
nt = 10000
ns = 10

nsave1 = 1
nwrte1 = 100

nsave2 = 10
nwrte2 = 100

h5file = "test.hdf5"


tab = ImplicitMidpoint()
ode = odeensemble()


# @testset "$(rpad("Serial Simulation",80))" begin

    # sim1 = Simulation(ode, tab, Δt, "Harmonic Oscillator Test 1 (Serial)", h5file, nt; nsave=nsave1, nwrite=nwrte1)
    # run!(sim1)
    # @test isfile(h5file)
    # ssol1 = SolutionODE(h5file)
    # @test ssol1.ntime == nt
    # @test ssol1.nsave == nsave1
    # rm(h5file)

    # sim2 = Simulation(ode, tab, Δt, "Harmonic Oscillator Test 2 (Serial)", h5file, nt; nsave=nsave2, nwrite=nwrte2)
    # run!(sim2)
    # @test isfile(h5file)
    # ssol2 = SolutionODE(h5file)
    # @test ssol2.ntime == nt
    # @test ssol2.nsave == nsave2
    # rm(h5file)

# end


# @testset "$(rpad("Parallel Simulation (" * string(Threads.nthreads()) * " Threads)",80))" begin

    # sim1 = ParallelSimulation(ode, tab, Δt, "Harmonic Oscillator Test 1 (" * string(Threads.nthreads()) * " Threads)", h5file, nt; nsave=nsave1, nwrite=nwrte1)
    # run!(sim1)
    # @test isfile(h5file)
    # psol1 = SolutionODE(h5file)
    # @test psol1.ntime == nt
    # @test psol1.nsave == nsave1
    # rm(h5file)

    # sim2 = ParallelSimulation(ode, tab, Δt, "Harmonic Oscillator Test 2 (" * string(Threads.nthreads()) * " Threads)", h5file, nt; nsave=nsave2, nwrite=nwrte2)
    # run!(sim2)
    # @test isfile(h5file)
    # psol2 = SolutionODE(h5file)
    # @test psol2.ntime == nt
    # @test psol2.nsave == nsave2
    # rm(h5file)

# end


# @test psol1.q == ssol1.q
# @test psol2.q == ssol2.q

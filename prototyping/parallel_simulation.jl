
using GeometricIntegrators
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using Test

set_config(:nls_atol_break, Inf)
set_config(:nls_rtol_break, Inf)
set_config(:nls_stol_break, Inf)


Δt = 0.1
nt = 100000
ns = 100

nsave = 1
nwrte = 100000

h5file = "test.hdf5"

ode = harmonic_oscillator_ode(vcat(rand(1,ns), zeros(1,ns)))
#tab = getTableauImplicitMidpoint()
tab = getTableauERK4()


sim = Simulation(ode, tab, Δt, "Serial Harmonic Oscillator Test", h5file, nt; nsave=nsave, nwrite=nwrte)
run!(sim)
rm(h5file)

sim = Simulation(ode, tab, Δt, "Serial Harmonic Oscillator Test", h5file, nt; nsave=nsave, nwrite=nwrte)
@time run!(sim)
rm(h5file)

sim = ParallelSimulation(ode, tab, Δt, "Parallel Harmonic Oscillator Test", h5file, nt; nsave=nsave, nwrite=nwrte)
run!(sim)
rm(h5file)

sim = ParallelSimulation(ode, tab, Δt, "Parallel Harmonic Oscillator Test", h5file, nt; nsave=nsave, nwrite=nwrte)
@time run!(sim)
rm(h5file)

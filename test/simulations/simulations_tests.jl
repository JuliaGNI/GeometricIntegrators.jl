
using GeometricIntegrators
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using Test

set_config(:nls_atol_break, Inf)
set_config(:nls_rtol_break, Inf)
set_config(:nls_stol_break, Inf)


const Δt = 0.1
const nt = 100000

const nsave1 = 1
const nwrte1 = 1000

const nsave2 = 10
const nwrte2 = 1000

const h5file = "test.hdf5"


ode = harmonic_oscillator_ode()
tab = getTableauImplicitMidpoint()



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

using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using SimpleSolvers
using Test

SimpleSolvers.set_config(:nls_stol_break, Inf)

using GeometricProblems.HarmonicOscillator: reference_solution

sode = harmonic_oscillator_sode()


sint = Integrator(sode, LieA())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

sintc = IntegratorComposition(sode, LieA())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, LieB())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

sintc = IntegratorComposition(sode, LieB())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, Strang())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

sintc = IntegratorComposition(sode, Strang())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, StrangA())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

sintc = IntegratorComposition(sode, StrangA())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, StrangB())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

sintc = IntegratorComposition(sode, StrangB())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint1 = Integrator(sode, Strang())
ssol1 = integrate(sode, sint)
sint2 = Integrator(sode, StrangA())
ssol2 = integrate(sode, sint)
sint3 = Integrator(sode, StrangB())
ssol3 = integrate(sode, sint)
@test ssol1.q == ssol2.q
@test all(ssol1.q .≈ ssol3.q)

sint = Integrator(sode, McLachlan2())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-4

sintc = IntegratorComposition(sode, McLachlan2())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, McLachlan4())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-4

sintc = IntegratorComposition(sode, McLachlan4())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, TripleJump())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-6

sintc = IntegratorComposition(sode, TripleJump())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q

sint = Integrator(sode, SuzukiFractal())
ssol = integrate(sode, sint)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-7

sintc = IntegratorComposition(sode, SuzukiFractal())
ssolc = integrate(sode, sintc)
@test ssol.q == ssolc.q


DT = datatype(sode)
D  = length(sode.ics.q)

ints_glrk1 = (IntegratorConstructor(DT, D, TableauGauss(1)), IntegratorConstructor(DT, D, TableauGauss(1)))
ints_erk4  = (IntegratorConstructor(DT, D, TableauRK4()), IntegratorConstructor(DT, D, TableauRK4()))

sint = IntegratorComposition(sode, ints_erk4, LieA())
ssol = integrate(sode, sint)
# println(relative_maximum_error(ssol.q, reference_solution))
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

sint = IntegratorComposition(sode, ints_glrk1, Strang())
ssol = integrate(sode, sint)
# println(relative_maximum_error(ssol.q, reference_solution))
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

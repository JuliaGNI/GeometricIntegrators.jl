
using GeometricBase.Config
using GeometricBase.Utils
using GeometricEquations
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricProblems.HarmonicOscillator
using SimpleSolvers
using Test

SimpleSolvers.set_config(:nls_stol_break, Inf)

using GeometricProblems.HarmonicOscillator: reference_solution, nt

sode = harmonic_oscillator_sode()


sint = Integrator(sode, TableauLieA())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

sintc = IntegratorComposition(sode, TableauLieA())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauLieB())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

sintc = IntegratorComposition(sode, TableauLieB())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauStrang())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

sintc = IntegratorComposition(sode, TableauStrang())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauStrangA())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

sintc = IntegratorComposition(sode, TableauStrangA())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauStrangB())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

sintc = IntegratorComposition(sode, TableauStrangB())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint1 = Integrator(sode, TableauStrang())
ssol1 = integrate(sode, sint, nt)
sint2 = Integrator(sode, TableauStrangA())
ssol2 = integrate(sode, sint, nt)
sint3 = Integrator(sode, TableauStrangB())
ssol3 = integrate(sode, sint, nt)
@test ssol1.q == ssol2.q
@test ssol1.q ≈  ssol3.q

sint = Integrator(sode, TableauMcLachlan2())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 1E-4

sintc = IntegratorComposition(sode, TableauMcLachlan2())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauMcLachlan4())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-4

sintc = IntegratorComposition(sode, TableauMcLachlan4())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauTripleJump())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-6

sintc = IntegratorComposition(sode, TableauTripleJump())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauSuzukiFractal())
ssol = integrate(sode, sint, nt)
@test relative_maximum_error(ssol.q, reference_solution) < 5E-7

sintc = IntegratorComposition(sode, TableauSuzukiFractal())
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q


DT = datatype(sode)
D  = length(sode.ics.q)

ints_glrk1 = (IntegratorConstructor(DT, D, TableauGauss(1)), IntegratorConstructor(DT, D, TableauGauss(1)))
ints_erk4  = (IntegratorConstructor(DT, D, TableauRK4()), IntegratorConstructor(DT, D, TableauRK4()))

sint = IntegratorComposition(sode, ints_erk4, TableauLieA())
ssol = integrate(sode, sint, nt)
# println(relative_maximum_error(ssol.q, reference_solution))
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

sint = IntegratorComposition(sode, ints_glrk1, TableauStrang())
ssol = integrate(sode, sint, nt)
# println(relative_maximum_error(ssol.q, reference_solution))
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

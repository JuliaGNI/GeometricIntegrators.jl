
set_config(:nls_solver, NewtonSolver)
set_config(:jacobian_autodiff, false)

function rel_err(sol, ref)
    maximum(abs.((sol.d[:,end] .- ref) ./ ref))
end

Δt  = 0.1
nt  = 10

refq = A * sin(ω * Δt * nt + ϕ)
refp = ω * Δt * nt * A * cos(ω * Δt * nt + ϕ)
refx = [refq, refp]

ode  = oscillator_ode()
pode = oscillator_pode()
sode = oscillator_sode()
iode = oscillator_iode()
idae = oscillator_idae()
pdae = oscillator_pdae()


@test typeof(Integrator(ode, getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
@test typeof(Integrator(ode, getTableauCrouzeix(), Δt)) <: IntegratorDIRK
@test typeof(Integrator(ode, getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK


### ERK Integrators ###

int = Integrator(ode, getTableauExplicitEuler(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-2

int = Integrator(ode, getTableauExplicitMidpoint(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 1E-3

int = Integrator(ode, getTableauERK4(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-7

### IRK Integrators ###

int = Integrator(ode, getTableauCrouzeix(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-5

int = Integrator(ode, getTableauImplicitEuler(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-2

int = Integrator(ode, getTableauImplicitMidpoint(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-4

int = IntegratorFIRK(ode, getTableauGLRK(1), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-4

int = IntegratorFIRK(ode, getTableauGLRK(2), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 1E-7

int = IntegratorFIRK(ode, getTableauGLRK(3), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 1E-11

int = IntegratorFIRK(ode, getTableauGLRK(4), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 1E-15

int = IntegratorFIRK(ode, getTableauGLRK(5), Δt)

@test rel_err(sol.q, refx) < 1E-15

int = IntegratorFIRK(ode, getTableauGLRK(6), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 1E-15

int = IntegratorFIRK(ode, getTableauGLRK(7), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 5E-15

int = IntegratorFIRK(ode, getTableauSRK3(), Δt)
sol = integrate(int, nt)

@test rel_err(sol.q, refx) < 1E-7

### PRK Integrators ###

pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
psol = integrate(pint, nt)

@test rel_err(psol.q, refq) < 5E-2
@test rel_err(psol.p, refp) < 1E-3

pint = Integrator(pode, getTableauSymplecticEulerB(), Δt)
psol = integrate(pint, nt)

@test rel_err(psol.q, refq) < 5E-2
@test rel_err(psol.p, refp) < 1E-3

pint = Integrator(pode, TableauEPRK(:prk4, 4, getTableauERK4().q), Δt)
psol = integrate(pint, nt)

@test rel_err(psol.q, refq) < 5E-7
@test rel_err(psol.p, refp) < 5E-7

pint = Integrator(pode, TableauIPRK(:pglrk, 2, getCoefficientsGLRK(1)), Δt)
psol = integrate(pint, nt)

@test rel_err(psol.q, refq) < 5E-4
@test rel_err(psol.p, refp) < 5E-4

pint = Integrator(pode, TableauIPRK(:pglrk, 4, getCoefficientsGLRK(2)), Δt)
psol = integrate(pint, nt)

@test rel_err(psol.q, refq) < 1E-7
@test rel_err(psol.p, refp) < 1E-7

pint = Integrator(pode, TableauIPRK(:pglrk, 6, getCoefficientsGLRK(3)), Δt)
psol = integrate(pint, nt)

@test rel_err(psol.q, refq) < 1E-11
@test rel_err(psol.p, refp) < 1E-11

### VPRK Integrators ###

iint = Integrator(iode, TableauVPRK(:pglrk, 2, getCoefficientsGLRK(1), -1), Δt)
isol = integrate(iint, nt)

@test rel_err(isol.q, refx) < 5E-4

iint = Integrator(iode, TableauVPRK(:pglrk, 4, getCoefficientsGLRK(2), +1), Δt)
isol = integrate(iint, nt)

@test rel_err(isol.q, refx) < 1E-7

iint = Integrator(iode, TableauVPRK(:pglrk, 6, getCoefficientsGLRK(3), -1), Δt)
isol = integrate(iint, nt)

@test rel_err(isol.q, refx) < 1E-11

vint = Integrator(iode, getTableauVPLobIIIA2(), Δt)
vsol = integrate(vint, nt)

@test rel_err(vsol.q, refx) < 5E-3

vint = Integrator(iode, getTableauVPLobIIIA3(), Δt)
vsol = integrate(vint, nt)

@test rel_err(vsol.q, refx) < 5E-4

vint = Integrator(iode, getTableauVPLobIIIA4(), Δt)
vsol = integrate(vint, nt)

@test rel_err(vsol.q, refx) < 1E-7

vint = Integrator(iode, getTableauVPLobIIIB2(), Δt)
vsol = integrate(vint, nt)

@test rel_err(vsol.q, refx) < 5E-3

vint = Integrator(iode, getTableauVPLobIIIB3(), Δt)
vsol = integrate(vint, nt)

@test rel_err(vsol.q, refx) < 5E-4

vint = Integrator(iode, getTableauVPLobIIIB4(), Δt)
vsol = integrate(vint, nt)

@test rel_err(vsol.q, refx) < 1E-7

vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 5E-4

vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(2), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-7

vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(3), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-11

vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 5E-4

vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(2), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-7

vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(3), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-11

vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 5E-4

vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(2), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-7

vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(3), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-11

vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 5E-4

vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(2), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-7

vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(3), Δt)
isol = integrate(vint, nt)

@test rel_err(isol.q, refx) < 1E-11

vintV1 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(1), Δt)
isolV1 = integrate(vintV1, nt)

@test rel_err(isolV1.q, refx) < 5E-4

vintV2 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(2), Δt)
isolV2 = integrate(vintV2, nt)

@test rel_err(isolV2.q, refx) < 1E-7

vintV3 = IntegratorVPRKpVariational(iode, getTableauVPGLRK(3), Δt)
isolV3 = integrate(vintV3, nt)

@test rel_err(isolV3.q, refx) < 1E-11

vintQ1 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(1), Δt)
isolQ1 = integrate(vintQ1, nt)

@test rel_err(isolQ1.q, refx) < 5E-4

vintQ2 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(2), Δt)
isolQ2 = integrate(vintQ2, nt)

@test rel_err(isolQ2.q, refx) < 1E-7

vintQ3 = IntegratorVPRKpVariationalQ(iode, getTableauVPGLRK(3), Δt)
isolQ3 = integrate(vintQ3, nt)

@test rel_err(isolQ3.q, refx) < 1E-11

vintP1 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(1), Δt)
isolP1 = integrate(vintP1, nt)

@test rel_err(isolP1.q, refx) < 5E-4

vintP2 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(2), Δt)
isolP2 = integrate(vintP2, nt)

@test rel_err(isolP2.q, refx) < 1E-7

vintP3 = IntegratorVPRKpVariationalP(iode, getTableauVPGLRK(3), Δt)
isolP3 = integrate(vintP3, nt)

@test rel_err(isolP3.q, refx) < 1E-11

@test rel_err(isolV1.q, isolP1.q[:,end]) < 1E-15
@test rel_err(isolV2.q, isolP2.q[:,end]) < 1E-15
@test rel_err(isolV3.q, isolP3.q[:,end]) < 1E-15


### Special Integrators ###

pgint = IntegratorPGLRK(iode, getCoefficientsPGLRK(2), Δt)
pgsol = integrate(pgint, nt)

@test rel_err(pgsol.q, refx) < 1E-5


### VPARK Integrator ###

dint = Integrator(idae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-4

dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-8

dint = Integrator(idae, getTableauGLRKpSymplectic(1), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-4

dint = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-8

dint = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 8E-4

dint = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 2E-4


### VSPARK Integrators ###

dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-4

dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-8

dint = Integrator(idae, getTableauGLRKpSymmetric(1), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-4

dint = Integrator(idae, getTableauGLRKpSymmetric(2), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-8

dint = Integrator(idae, getTableauLobIIIAIIIB2pSymmetric(), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 8E-3



### HPARK Integrator ###

dint = Integrator(pdae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-4

dint = Integrator(pdae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-8

dint = Integrator(pdae, getTableauGLRKpSymplectic(1), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-4

dint = Integrator(pdae, getTableauGLRKpSymplectic(2), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 5E-8

dint = Integrator(pdae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 2E-2

dint = Integrator(pdae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
dsol = integrate(dint, nt)

@test rel_err(dsol.q, refx) < 2E-4


### CGVI and DGVI Integrators ###

QGau4 = GaussLegendreQuadrature(4)
BGau4 = LagrangeBasis(nodes(QGau4))

cgint = IntegratorCGVI(iode, BGau4, QGau4, Δt)
cgsol = integrate(cgint, nt)

@test rel_err(cgsol.q, refx) < 1E-7

dgint = IntegratorDGVI(iode, BGau4, QGau4, Δt)
dgsol = integrate(dgint, nt)

@test rel_err(dgsol.q, refx) < 1E-7


### Splitting Integrators ###

sint = Integrator(sode, getTableauLieA(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 5E-2

sint = Integrator(sode, getTableauLieB(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 5E-2

sint = Integrator(sode, getTableauStrang(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 1E-3

sint = Integrator(sode, getTableauMcLachlan2(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 1E-4

sint = Integrator(sode, getTableauMcLachlan4(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 5E-4

sint = Integrator(sode, getTableauTripleJump(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 5E-6

sint = Integrator(sode, getTableauSuzukiFractal(), Δt)
ssol = integrate(sint, nt)

@test rel_err(ssol.q, refx) < 5E-7

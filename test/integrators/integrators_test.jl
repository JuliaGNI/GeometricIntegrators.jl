
set_config(:nls_solver, NewtonSolver)

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

# vint = Integrator(iode, getCoefficientsPGLRK(2), Δt)
# vsol = integrate(vint, nt)

# TODO Add PDAE/PARK test.

# dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2p, glrk2.q, glrk2.q), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2p, glrk2.q, glrk2.q), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauGLRKpSymmetric(2), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauLobIIIAIIIB2pSymmetric(), Δt)
# dsol = integrate(dint, nt)

# dint = Integrator(idae, getTableauLobIIIAIIIB3pSymmetric(), Δt)
# dsol = integrate(dint, nt)


### Special Integrators ###

pgint = IntegratorPGLRK(iode, getCoefficientsPGLRK(2), Δt)
pgsol = integrate(pgint, nt)

@test rel_err(pgsol.q, refx) < 1E-5


### CGVI and DGVI Integrators ###

QGau4 = GaussLegendreQuadrature(4)
BGau4 = LagrangeBasis(nodes(QGau4))
cgint = IntegratorCGVI(iode, BGau4, QGau4, Δt)
cgsol = integrate(cgint, nt)

@test rel_err(cgsol.q, refx) < 1E-7


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

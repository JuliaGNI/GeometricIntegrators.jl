
set_config(:nls_solver, NewtonSolver)
set_config(:jacobian_autodiff, false)

iode = oscillator_iode()


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

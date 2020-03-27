
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem
using Test

using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: Δt

int = get_config(:ig_interpolation)

ode  = lotka_volterra_2d_ode()
pode = lotka_volterra_2d_pode()
iode = lotka_volterra_2d_iode()
vode = lotka_volterra_2d_vode()

pdae = lotka_volterra_2d_pdae()
idae = lotka_volterra_2d_idae()
vdae = lotka_volterra_2d_vdae()


@test InitialGuessODE(int,  ode, Δt) == InitialGuessODE{eltype( ode.q₀), ndims( ode)}(int,  ode.v, Δt)
@test InitialGuessODE(int, iode, Δt) == InitialGuessODE{eltype(iode.q₀), ndims(iode)}(int, iode.v, Δt)
@test InitialGuessODE(int, vode, Δt) == InitialGuessODE{eltype(vode.q₀), ndims(vode)}(int, vode.v, Δt)

@test InitialGuessPODE(int, pode, Δt) == InitialGuessPODE{eltype(pode.q₀), ndims(pode)}(int, pode.v, pode.f, Δt)
@test InitialGuessPODE(int, iode, Δt) == InitialGuessPODE{eltype(iode.q₀), ndims(iode)}(int, iode.v, iode.f, Δt)
@test InitialGuessPODE(int, vode, Δt) == InitialGuessPODE{eltype(vode.q₀), ndims(vode)}(int, vode.v, vode.f, Δt)

@test InitialGuessPODE(int, pdae, Δt) == InitialGuessPODE{eltype(pdae.q₀), ndims(pdae)}(int, pdae.v, pdae.f, Δt)
@test InitialGuessPODE(int, idae, Δt) == InitialGuessPODE{eltype(idae.q₀), ndims(idae)}(int, idae.v, idae.f, Δt)
@test InitialGuessPODE(int, vdae, Δt) == InitialGuessPODE{eltype(vdae.q₀), ndims(vdae)}(int, vdae.v, vdae.f, Δt)

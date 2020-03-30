
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


@test InitialGuess(int,  ode, Δt) == InitialGuessODE{eltype(ode.q₀), ndims(ode)}(int,  ode.v, Δt)

@test InitialGuess(int, iode, Δt) == InitialGuessIODE{eltype(iode.q₀), ndims(iode)}(int, iode.v, iode.f, Δt)
@test InitialGuess(int, idae, Δt) == InitialGuessIODE{eltype(idae.q₀), ndims(idae)}(int, idae.v, idae.f, Δt)

@test InitialGuess(int, vode, Δt) == InitialGuessIODE{eltype(vode.q₀), ndims(vode)}(int, vode.v, vode.f, Δt)
@test InitialGuess(int, vdae, Δt) == InitialGuessIODE{eltype(vdae.q₀), ndims(vdae)}(int, vdae.v, vdae.f, Δt)

@test InitialGuess(int, pode, Δt) == InitialGuessPODE{eltype(pode.q₀), ndims(pode)}(int, pode.v, pode.f, Δt)
@test InitialGuess(int, pdae, Δt) == InitialGuessPODE{eltype(pdae.q₀), ndims(pdae)}(int, pdae.v, pdae.f, Δt)

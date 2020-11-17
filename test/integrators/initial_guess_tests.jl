
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
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
@test InitialGuessODE(int, iode, Δt) == InitialGuessODE{eltype(iode.q₀), ndims(iode)}(int, iode.v̄, Δt)
@test InitialGuessODE(int, vode, Δt) == InitialGuessODE{eltype(vode.q₀), ndims(vode)}(int, vode.v̄, Δt)

@test InitialGuess(int, ode, Δt) == InitialGuessODE{eltype(ode.q₀), ndims(ode)}(int, ode.v, Δt)

@test InitialGuess(int, iode, Δt) == InitialGuessIODE{eltype(iode.q₀), ndims(iode)}(int, iode.v̄, iode.f̄, Δt)
@test InitialGuess(int, idae, Δt) == InitialGuessIODE{eltype(idae.q₀), ndims(idae)}(int, idae.v̄, idae.f̄, Δt)

@test InitialGuess(int, vode, Δt) == InitialGuessIODE{eltype(vode.q₀), ndims(vode)}(int, vode.v̄, vode.f̄, Δt)
@test InitialGuess(int, vdae, Δt) == InitialGuessIODE{eltype(vdae.q₀), ndims(vdae)}(int, vdae.v̄, vdae.f̄, Δt)

@test InitialGuess(int, pode, Δt) == InitialGuessPODE{eltype(pode.q₀), ndims(pode)}(int, pode.v, pode.f, Δt)
@test InitialGuess(int, pdae, Δt) == InitialGuessPODE{eltype(pdae.q₀), ndims(pdae)}(int, pdae.v̄, pdae.f̄, Δt)


# Reference Solution

ref_prev = integrate(ode, getTableauERK4(), -Δt/10, 10)
ref_next = integrate(ode, getTableauERK4(), +Δt/10, 10)

tₚ = ref_prev.t[end]
qₚ = ref_prev.q[:,end]
vₚ = zero(qₚ)

tₙ = ref_next.t[end]
qₙ = ref_next.q[:,end]
vₙ = zero(qₙ)

ode.v(tₚ, qₚ, vₚ, ode.parameters)
ode.v(tₙ, qₙ, vₙ, ode.parameters)


# InitialGuessODE

igode = InitialGuessODE{eltype(ode.q₀), ndims(ode)}(int, (t,q,v) -> ode.v(t, q, v, ode.parameters), Δt)

t₀ = ode.t₀
q₀ = ode.q₀
v₀ = zero(q₀)

t₁ = ode.t₀ - Δt
q₁ = zero(q₀)
v₁ = zero(v₀)

t₂ = one(Δt)
q₂ = zero(q₀)
v₂ = zero(v₀)

ode.v(t₀, q₀, v₀, ode.parameters)

initialize!(igode, t₀, q₀, v₀, t₁, q₁, v₁)
evaluate!(igode, q₁, v₁, q₀, v₀, q₂, v₂, t₂)

@test q₁ ≈ qₚ atol=1E-5
@test v₁ ≈ vₚ atol=1E-5

@test q₂ ≈ qₙ atol=1E-4
@test v₂ ≈ vₙ atol=1E-2


# InitialGuessIODE

igiode = InitialGuessIODE{eltype(iode.q₀), ndims(iode)}(int, (t,q,v) -> iode.v̄(t, q, v, iode.parameters), (t,q,p,v) -> iode.f̄(t, q, p, v, iode.parameters), Δt)

t₀ = iode.t₀
q₀ = iode.q₀
p₀ = iode.p₀
v₀ = zero(q₀)
f₀ = zero(p₀)

t₁ = iode.t₀ - Δt
q₁ = zero(q₀)
p₁ = zero(p₀)
v₁ = zero(v₀)
f₁ = zero(f₀)

t₂ = one(Δt)
q₂ = zero(q₀)
p₂ = zero(p₀)
v₂ = zero(v₀)
f₂ = zero(f₀)

iode.v̄(t₀, q₀, v₀, iode.parameters)

initialize!(igiode, t₀, q₀, p₀, v₀, f₀, t₁, q₁, p₁, v₁, f₁)
evaluate!(igiode, q₁, p₁, v₁, f₁, q₀, p₀, v₀, f₀, q₂, v₂, t₂)

@test q₁ ≈ qₚ atol=1E-5
@test v₁ ≈ vₚ atol=1E-5

@test q₂ ≈ qₙ atol=1E-4
@test v₂ ≈ vₙ atol=1E-2


# InitialGuessPODE

igpode = InitialGuessPODE{eltype(pode.q₀), ndims(pode)}(int,
            (t,q,p,v) -> pode.v(t, q, p, v, pode.parameters),
            (t,q,p,v) -> pode.f(t, q, p, v, pode.parameters), Δt)

t₀ = pode.t₀
q₀ = pode.q₀
p₀ = pode.p₀
v₀ = zero(q₀)
f₀ = zero(p₀)

t₁ = pode.t₀ - Δt
q₁ = zero(q₀)
p₁ = zero(p₀)
v₁ = zero(v₀)
f₁ = zero(f₀)

t₂ = one(Δt)
q₂ = zero(q₀)
p₂ = zero(p₀)
v₂ = zero(v₀)
f₂ = zero(f₀)

pode.v(t₀, q₀, p₀, v₀, pode.parameters)
pode.f(t₀, q₀, p₀, f₀, pode.parameters)

initialize!(igpode, t₀, q₀, p₀, v₀, f₀, t₁, q₁, p₁, v₁, f₁)
evaluate!(igpode, q₁, p₁, v₁, f₁, q₀, p₀, v₀, f₀, q₂, v₂, t₂)

@test q₁ ≈ qₚ atol=1E-5
@test v₁ ≈ vₚ atol=1E-5

@test q₂ ≈ qₙ atol=1E-4
@test v₂ ≈ vₙ atol=1E-2

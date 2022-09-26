using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using GeometricProblems.LotkaVolterra2d: Δt
using Test

using GeometricEquations: _get_v, _get_f, _get_v̄, _get_f̄

int = get_config(:ig_extrapolation)

const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

ode  = lotka_volterra_2d_ode(q₀; parameters=params)
pode = lotka_volterra_2d_pode(q₀; parameters=params)
hode = lotka_volterra_2d_hode(q₀; parameters=params)
iode = lotka_volterra_2d_iode(q₀; parameters=params)
lode = lotka_volterra_2d_lode(q₀; parameters=params)

dae  = lotka_volterra_2d_dae(q₀; parameters=params)
pdae = lotka_volterra_2d_pdae(q₀; parameters=params)
hdae = lotka_volterra_2d_hdae(q₀; parameters=params)
idae = lotka_volterra_2d_idae(q₀; parameters=params)
ldae = lotka_volterra_2d_ldae(q₀; parameters=params)


@test InitialGuessODE(int,  ode) == InitialGuessODE(int(0.0, Δt), _get_v̄(equation( ode), params), Δt)
@test InitialGuessODE(int, iode) == InitialGuessODE(int(0.0, Δt), _get_v̄(equation(iode), params), Δt)
@test InitialGuessODE(int, lode) == InitialGuessODE(int(0.0, Δt), _get_v̄(equation(lode), params), Δt)

@test InitialGuess(int, ode) == InitialGuessODE(int(0.0, Δt), _get_v̄(equation(ode), params), Δt)
@test InitialGuess(int, dae) == InitialGuessODE(int(0.0, Δt), _get_v̄(equation(dae), params), Δt)

@test InitialGuess(int, hode) == InitialGuessPODE(int(0.0, Δt), _get_v̄(equation(hode), params), _get_f̄(equation(hode), params), Δt)
@test InitialGuess(int, hdae) == InitialGuessPODE(int(0.0, Δt), _get_v̄(equation(hdae), params), _get_f̄(equation(hdae), params), Δt)

@test InitialGuess(int, iode) == InitialGuessIODE(int(0.0, Δt), _get_v̄(equation(iode), params), _get_f̄(equation(iode), params), Δt)
@test InitialGuess(int, idae) == InitialGuessIODE(int(0.0, Δt), _get_v̄(equation(idae), params), _get_f̄(equation(idae), params), Δt)

@test InitialGuess(int, pode) == InitialGuessPODE(int(0.0, Δt), _get_v̄(equation(pode), params), _get_f̄(equation(pode), params), Δt)
@test InitialGuess(int, pdae) == InitialGuessPODE(int(0.0, Δt), _get_v̄(equation(pdae), params), _get_f̄(equation(pdae), params), Δt)

@test InitialGuess(int, lode) == InitialGuessIODE(int(0.0, Δt), _get_v̄(equation(lode), params), _get_f̄(equation(lode), params), Δt)
@test InitialGuess(int, ldae) == InitialGuessIODE(int(0.0, Δt), _get_v̄(equation(ldae), params), _get_f̄(equation(ldae), params), Δt)


# Reference Solution

ref_prev = integrate(similar(ode; tspan=(tspan(ode)[begin], tspan(ode)[begin]-tstep(ode)), tstep=-tstep(ode)), TableauGauss(8))
ref_next = integrate(similar(ode; tspan=(tspan(ode)[begin], tspan(ode)[begin]+tstep(ode)), tstep=+tstep(ode)), TableauGauss(8))

tₚ = ref_prev.t[end]
qₚ = ref_prev.q[end]
vₚ = zero(qₚ)

tₙ = ref_next.t[end]
qₙ = ref_next.q[end]
vₙ = zero(qₙ)

equation(ode).v(tₚ, qₚ, vₚ, parameters(ode))
equation(ode).v(tₙ, qₙ, vₙ, parameters(ode))


# InitialGuessODE

igode = InitialGuessODE(int, _get_v(equation(ode), parameters(ode)), Δt)

t₀ = tbegin(ode)
q₀ = ode.ics.q
v₀ = zero(q₀)

t₁ = tbegin(ode) - Δt
q₁ = zero(q₀)
v₁ = zero(v₀)

t₂ = one(Δt)
q₂ = zero(q₀)
v₂ = zero(v₀)

equation(ode).v(t₀, q₀, v₀, parameters(ode))

initialize!(igode, t₀, q₀, v₀, t₁, q₁, v₁)
evaluate!(igode, q₁, v₁, q₀, v₀, q₂, v₂, t₂)

# println("IG-ODE")
# println(q₁ .- qₚ)
# println(v₁ .- vₚ)
# println(q₂ .- qₙ)
# println(v₂ .- vₙ)
# println()

@test q₁ ≈ qₚ atol=1E-14
@test v₁ ≈ vₚ atol=1E-14

@test q₂ ≈ qₙ atol=1E-8
@test v₂ ≈ vₙ atol=1E-5


# InitialGuessIODE

igiode = InitialGuessIODE(int, _get_v̄(equation(iode), parameters(ode)), _get_f̄(equation(iode), parameters(ode)), Δt)

t₀ = tbegin(iode)
q₀ = iode.ics.q
p₀ = iode.ics.p
v₀ = zero(q₀)
f₀ = zero(p₀)

t₁ = tbegin(iode) - Δt
q₁ = zero(q₀)
p₁ = zero(p₀)
v₁ = zero(v₀)
f₁ = zero(f₀)

t₂ = one(Δt)
q₂ = zero(q₀)
p₂ = zero(p₀)
v₂ = zero(v₀)
f₂ = zero(f₀)

equation(iode).v̄(t₀, q₀, v₀, parameters(iode))

initialize!(igiode, t₀, q₀, p₀, v₀, f₀, t₁, q₁, p₁, v₁, f₁)
evaluate!(igiode, q₁, p₁, v₁, f₁, q₀, p₀, v₀, f₀, q₂, v₂, t₂)

# println("IG-IODE")
# println(q₁ .- qₚ)
# println(v₁ .- vₚ)
# println(q₂ .- qₙ)
# println(v₂ .- vₙ)
# println()

@test q₁ ≈ qₚ atol=1E-14
@test v₁ ≈ vₚ atol=1E-14

@test q₂ ≈ qₙ atol=1E-8
@test v₂ ≈ vₙ atol=1E-5


# InitialGuessPODE

igpode = InitialGuessPODE(int, _get_v(equation(pode), parameters(pode)), _get_f(equation(pode), parameters(pode)), Δt)

t₀ = tbegin(pode)
q₀ = pode.ics.q
p₀ = pode.ics.p
v₀ = zero(q₀)
f₀ = zero(p₀)

t₁ = tbegin(pode) - Δt
q₁ = zero(q₀)
p₁ = zero(p₀)
v₁ = zero(v₀)
f₁ = zero(f₀)

t₂ = one(Δt)
q₂ = zero(q₀)
p₂ = zero(p₀)
v₂ = zero(v₀)
f₂ = zero(f₀)

equation(pode).v(t₀, q₀, p₀, v₀, parameters(pode))
equation(pode).f(t₀, q₀, p₀, f₀, parameters(pode))

initialize!(igpode, t₀, q₀, p₀, v₀, f₀, t₁, q₁, p₁, v₁, f₁)
evaluate!(igpode, q₁, p₁, v₁, f₁, q₀, p₀, v₀, f₀, q₂, v₂, t₂)

# println("IG-PODE")
# println(q₁ .- qₚ)
# println(v₁ .- vₚ)
# println(q₂ .- qₙ)
# println(v₂ .- vₙ)
# println()

@test q₁ ≈ qₚ atol=1E-14
@test v₁ ≈ vₚ atol=1E-14

@test q₂ ≈ qₙ atol=1E-8
@test v₂ ≈ vₙ atol=1E-5

using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

ode  = odeproblem()
pode = podeproblem()
iode = degenerate_iodeproblem()


# Compute Reference Solution for ODEs

const Δt = timestep(ode)
const t₀ = initial_conditions(ode).t
const t₁ = t₀ + Δt
const t₂ = t₁ + Δt
const t₋ = t₀ - Δt
const tₚ = t₋
const tₙ = t₁
const tᵢ = tₙ

x₀ = initial_conditions(ode).q

k = parameters(ode).k
ω = parameters(ode).ω
A = sqrt(x₀[2]^2 / k + x₀[1]^2)
ϕ = asin(x₀[1] / A)

xₚ = exact_solution(t₀ - Δt, x₀, t₀, parameters(ode))
xₙ = exact_solution(t₀ + Δt, x₀, t₀, parameters(ode))


# Create ODE Solution Arrays

x₁ = zero(x₀)
x₂ = zero(x₀)
xᵢ = zero(x₀)

ẋₚ = zero(x₀)
ẋ₀ = zero(x₀)
ẋ₁ = zero(x₀)
ẋ₂ = zero(x₀)
ẋₙ = zero(x₀)
ẋᵢ = zero(x₀)

functions(ode).v(ẋₚ, tₚ, xₚ, parameters(ode))
functions(ode).v(ẋ₀, t₀, x₀, parameters(ode))
functions(ode).v(ẋₙ, tₙ, xₙ, parameters(ode))


# Hermite Extrapolation

extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, tᵢ, xᵢ, ẋᵢ, HermiteExtrapolation())

# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, ẋₙ, ẋᵢ .- ẋₙ)

@test xᵢ ≈ xₙ atol=1E-5
@test ẋᵢ ≈ ẋₙ atol=1E-4

extrapolate!(tₚ, xₚ, t₀, x₀, t₁, x₁, ode, HermiteExtrapolation())

@test x₁ == xᵢ

@test extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, t₁, x₁, HermiteExtrapolation()) == xᵢ
@test extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, HermiteExtrapolation()) == (xᵢ, ẋᵢ)
@test extrapolate!(tₚ, xₚ, t₀, x₀, t₁, x₁, functions(ode).v, parameters(ode), HermiteExtrapolation()) == xᵢ
@test extrapolate!(tₚ, xₚ, t₀, x₀, t₁, x₁, ẋ₁, functions(ode).v, parameters(ode), HermiteExtrapolation()) == (xᵢ, ẋᵢ)
@test extrapolate!(tₚ, xₚ, t₀, x₀, t₁, x₁, ode, HermiteExtrapolation()) == xᵢ
@test extrapolate!(tₚ, xₚ, t₀, x₀, t₁, x₁, ẋ₁, ode, HermiteExtrapolation()) == (xᵢ, ẋᵢ)


# solution and history tuples
sol = (t = tᵢ, q = xᵢ)
history = (t = [t₀], q = [x₀])

# Euler Extrapolation for ODEs

extrapolate!(sol, history, ode, EulerExtrapolation(0))
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test sol.q ≈ xₙ atol=1E-1

extrapolate!(sol, history, ode, EulerExtrapolation(1))
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test sol.q ≈ xₙ atol=1E-2

extrapolate!(sol, history, ode, EulerExtrapolation(2))
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test sol.q ≈ xₙ atol=1E-4

extrapolate!(sol, history, ode, EulerExtrapolation(3))
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test sol.q ≈ xₙ atol=1E-6

extrapolate!(sol, history, ode, EulerExtrapolation(4))
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test sol.q ≈ xₙ atol=1E-8

extrapolate!(sol, history, ode, EulerExtrapolation(5))
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test sol.q ≈ xₙ atol=1E-10


# Midpoint Extrapolation for ODEs

extrapolate!(sol, history, ode, MidpointExtrapolation(0))
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test sol.q ≈ xₙ atol=1E-4

extrapolate!(sol, history, ode, MidpointExtrapolation(1))
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test sol.q ≈ xₙ atol=1E-8

extrapolate!(sol, history, ode, MidpointExtrapolation(2))
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test sol.q ≈ xₙ atol=1E-12

extrapolate!(sol, history, ode, MidpointExtrapolation(3))
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test sol.q ≈ xₙ atol=1E-15

extrapolate!(sol, history, ode, MidpointExtrapolation(4))
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test sol.q ≈ xₙ atol=1E-16

extrapolate!(sol, history, ode, MidpointExtrapolation(5))
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test sol.q ≈ xₙ atol=1E-15


# Create PODE Solution Arrays

q₀ = initial_conditions(pode).q
p₀ = initial_conditions(pode).p

qᵢ = zero(q₀)
pᵢ = zero(p₀)

q̇ₚ = zero(q₀)
q̇₀ = zero(q₀)
q̇ₙ = zero(q₀)
q̇ᵢ = zero(q₀)

ṗₚ = zero(p₀)
ṗ₀ = zero(p₀)
ṗₙ = zero(p₀)
ṗᵢ = zero(p₀)


# Compute Reference Solution for PODEs

qₚ = [xₚ[1]]
pₚ = [xₚ[2]]

qₙ = [xₙ[1]]
pₙ = [xₙ[2]]

functions(pode).v(q̇ₚ, tₚ, qₚ, pₚ, parameters(pode))
functions(pode).v(q̇₀, t₀, q₀, p₀, parameters(pode))
functions(pode).v(q̇ₙ, tₙ, qₙ, pₙ, parameters(pode))

functions(pode).f(ṗₚ, tₚ, qₚ, pₚ, parameters(pode))
functions(pode).f(ṗ₀, t₀, q₀, p₀, parameters(pode))
functions(pode).f(ṗₙ, tₙ, qₙ, pₙ, parameters(pode))

# solution and history tuples
sol = (t = tᵢ, q = qᵢ, p = pᵢ)
history = (t = [t₀], q = [q₀], p = [p₀])

# Midpoint Extrapolation for PODEs

extrapolate!(sol, history, pode, MidpointExtrapolation(0))
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-6
@test pᵢ ≈ pₙ atol=1E-4

extrapolate!(sol, history, pode, MidpointExtrapolation(1))
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-10
@test pᵢ ≈ pₙ atol=1E-8

extrapolate!(sol, history, pode, MidpointExtrapolation(2))
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-14
@test pᵢ ≈ pₙ atol=1E-12

extrapolate!(sol, history, pode, MidpointExtrapolation(3))
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

extrapolate!(sol, history, pode, MidpointExtrapolation(4))
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-16
@test pᵢ ≈ pₙ atol=1E-16

extrapolate!(sol, history, pode, MidpointExtrapolation(5))
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16


# Create IODE Solution Arrays

q₀ = initial_conditions(iode).q
p₀ = initial_conditions(iode).p

qᵢ = zero(q₀)
qₚ = zero(q₀)
qₙ = zero(q₀)

pᵢ = zero(p₀)
pₚ = zero(p₀)
pₙ = zero(p₀)

q̇ₚ = zero(q₀)
q̇₀ = zero(q₀)
q̇ₙ = zero(q₀)
q̇ᵢ = zero(q₀)

ṗₚ = zero(p₀)
ṗ₀ = zero(p₀)
ṗₙ = zero(p₀)
ṗᵢ = zero(p₀)


# Compute Reference Solution for IODEs

qₚ .= xₚ
qₙ .= xₙ

functions(iode).v̄(q̇ₚ, tₚ, qₚ, pₚ, parameters(iode))
functions(iode).v̄(q̇₀, t₀, q₀, p₀, parameters(iode))
functions(iode).v̄(q̇ₙ, tₙ, qₙ, pₙ, parameters(iode))

functions(iode).ϑ(pₚ, tₚ, qₚ, q̇ₚ, parameters(iode))
functions(iode).ϑ(pₙ, tₙ, qₙ, q̇ₙ, parameters(iode))

functions(iode).f̄(ṗₚ, tₚ, qₚ, q̇ₚ, parameters(iode))
functions(iode).f̄(ṗ₀, t₀, q₀, q̇₀, parameters(iode))
functions(iode).f̄(ṗₙ, tₙ, qₙ, q̇ₙ, parameters(iode))

# solution and history tuples
sol = (t = tᵢ, q = qᵢ, p = pᵢ)
history = (t = [t₀], q = [q₀], p = [p₀])

# Midpoint Extrapolation for IODEs

extrapolate!(sol, history, iode, MidpointExtrapolation(0))
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-4
@test pᵢ ≈ pₙ atol=1E-4

extrapolate!(sol, history, iode, MidpointExtrapolation(1))
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-8
@test pᵢ ≈ pₙ atol=1E-8

extrapolate!(sol, history, iode, MidpointExtrapolation(2))
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-12
@test pᵢ ≈ pₙ atol=1E-12

extrapolate!(sol, history, iode, MidpointExtrapolation(3))
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

extrapolate!(sol, history, iode, MidpointExtrapolation(4))
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-16
@test pᵢ ≈ pₙ atol=1E-16

extrapolate!(sol, history, iode, MidpointExtrapolation(5))
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

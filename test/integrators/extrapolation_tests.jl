
using GeometricIntegrators.Common
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricProblems.HarmonicOscillator
using Test

using GeometricIntegrators.Equations: _get_v
using GeometricIntegrators.Integrators: _euler_extrapolation_ode!, midpoint_extrapolation

const Δt = 0.1
const t₀ = 0.0
const t₁ = t₀ + Δt

ode  = harmonic_oscillator_ode()


# Create Solution Arrays

q₀ = ode.q₀[1]
qᵢ = zero(q₀)

vₚ = zero(q₀)
v₀ = zero(q₀)
vₙ = zero(q₀)
vᵢ = zero(q₀)


# Compute Reference Solution

ref_prev = integrate(ode, TableauGauss(8), -Δt, 1)
ref_next = integrate(ode, TableauGauss(8), +Δt, 1)

tₚ = ref_prev.t[end]
qₚ = ref_prev.q[end]

tₙ = ref_next.t[end]
qₙ = ref_next.q[end]


# Compute Vector Fields

ode.v(tₚ, qₚ, vₚ, ode.parameters)
ode.v(t₀, q₀, v₀, ode.parameters)
ode.v(tₙ, qₙ, vₙ, ode.parameters)


# Hermite Extrapolation

int = HermiteExtrapolation(zero(Δt), one(Δt), Δt)
evaluate!(int, qₚ, q₀, vₚ, v₀, 2one(Δt), qᵢ, vᵢ)

# println(qᵢ, qₙ, qᵢ .- qₙ)
# println(vᵢ, vₙ, vᵢ .- vₙ)

@test qᵢ ≈ qₙ atol=1E-5
@test vᵢ ≈ vₙ atol=1E-4


# Euler Extrapolation

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, 0)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-1

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, 1)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-2

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, 2)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-4

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, 3)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-6

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, 4)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-8

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, 5)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-10

for i in 0:5
    @test EulerExtrapolation(ode, i) == EulerExtrapolationODE(_get_v(ode), i)
    @test evaluate!(EulerExtrapolation(ode, i), t₀, t₁, q₀, qᵢ) == _euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, i)
    @test evaluate!(EulerExtrapolationODE(_get_v(ode), i), t₀, t₁, q₀, qᵢ) == _euler_extrapolation_ode!(_get_v(ode), t₀, t₁, q₀, qᵢ, i)
end


# Midpoint Extrapolation

midpoint_extrapolation(_get_v(ode), t₀, t₁, q₀::Vector, qᵢ::Vector, 0)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-4

midpoint_extrapolation(_get_v(ode), t₀, t₁, q₀::Vector, qᵢ::Vector, 1)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-8

midpoint_extrapolation(_get_v(ode), t₀, t₁, q₀::Vector, qᵢ::Vector, 2)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-12

midpoint_extrapolation(_get_v(ode), t₀, t₁, q₀::Vector, qᵢ::Vector, 3)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-16

midpoint_extrapolation(_get_v(ode), t₀, t₁, q₀::Vector, qᵢ::Vector, 4)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-16

midpoint_extrapolation(_get_v(ode), t₀, t₁, q₀::Vector, qᵢ::Vector, 5)
# println(qᵢ, qₙ, qᵢ .- qₙ)
@test qᵢ ≈ qₙ atol=1E-15


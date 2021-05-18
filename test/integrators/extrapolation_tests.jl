
using GeometricIntegrators.Common
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricProblems.HarmonicOscillator
using Test

using GeometricIntegrators.Equations: _get_v, _get_f, _get_v̄, _get_f̄
using GeometricIntegrators.Integrators: _euler_extrapolation_ode!, _midpoint_extrapolation_ode!,
                                        _midpoint_extrapolation_iode!, _midpoint_extrapolation_pode!

const Δt = 0.1
const t₀ = 0.0
const t₁ = t₀ + Δt

ode  = harmonic_oscillator_ode()
pode = harmonic_oscillator_pode()
iode = harmonic_oscillator_iode()


# Create ODE Solution Arrays

x₀ = ode.q₀[1]
xᵢ = zero(x₀)

ẋₚ = zero(x₀)
ẋ₀ = zero(x₀)
ẋₙ = zero(x₀)
ẋᵢ = zero(x₀)


# Compute Reference Solution for ODEs

ref_prev = integrate(ode, TableauGauss(8), -Δt, 1)
ref_next = integrate(ode, TableauGauss(8), +Δt, 1)

tₚ = ref_prev.t[end]
xₚ = ref_prev.q[end]

tₙ = ref_next.t[end]
xₙ = ref_next.q[end]

ode.v(tₚ, xₚ, ẋₚ, ode.parameters)
ode.v(t₀, x₀, ẋ₀, ode.parameters)
ode.v(tₙ, xₙ, ẋₙ, ode.parameters)


# Create PODE Solution Arrays

q₀ = pode.q₀[1]
p₀ = pode.p₀[1]

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

ref_prev = integrate(pode, PartitionedTableauGauss(8), -Δt, 1)
ref_next = integrate(pode, PartitionedTableauGauss(8), +Δt, 1)

qₚ = ref_prev.q[end]
pₚ = ref_prev.p[end]

qₙ = ref_next.q[end]
pₙ = ref_next.p[end]

pode.v(tₚ, qₚ, pₚ, q̇ₚ, ode.parameters)
pode.v(t₀, q₀, p₀, q̇₀, ode.parameters)
pode.v(tₙ, qₙ, pₙ, q̇ₙ, ode.parameters)

pode.f(tₚ, qₚ, pₚ, ṗₚ, ode.parameters)
pode.f(t₀, q₀, p₀, ṗ₀, ode.parameters)
pode.f(tₙ, qₙ, pₙ, ṗₙ, ode.parameters)



# Hermite Extrapolation

int = HermiteExtrapolation(zero(Δt), one(Δt), Δt)
evaluate!(int, xₚ, x₀, ẋₚ, ẋ₀, 2one(Δt), xᵢ, ẋᵢ)

# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, ẋₙ, ẋᵢ .- ẋₙ)

@test xᵢ ≈ xₙ atol=1E-5
@test ẋᵢ ≈ ẋₙ atol=1E-4


# Euler Extrapolation

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 0)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-1

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 1)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-2

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 2)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-4

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 3)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-6

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 4)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-8

_euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 5)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-10

for i in 0:5
    extrap1 = EulerExtrapolation(ode, i)
    extrap2 = EulerExtrapolationODE(_get_v(ode), i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, x₀, xᵢ) == _euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, x₀, xᵢ) == _euler_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, i)
end


# Midpoint Extrapolation for ODEs

_midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 0)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-4

_midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 1)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-8

_midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 2)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-12

_midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 3)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-16

_midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 4)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-16

_midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, 5)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-15

for i in 0:5
    extrap1 = MidpointExtrapolation(ode, i)
    extrap2 = MidpointExtrapolationODE(_get_v(ode), i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, x₀, xᵢ) == _midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, x₀, xᵢ) == _midpoint_extrapolation_ode!(_get_v(ode), t₀, t₁, x₀, xᵢ, i)
end


# Midpoint Extrapolation for PODEs

_midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, 0)
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-4
@test pᵢ ≈ pₙ atol=1E-4

_midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, 1)
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-8
@test pᵢ ≈ pₙ atol=1E-8

_midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, 2)
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-12
@test pᵢ ≈ pₙ atol=1E-12

_midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, 3)
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-16
@test pᵢ ≈ pₙ atol=1E-16

_midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, 4)
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-16
@test pᵢ ≈ pₙ atol=1E-16

_midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, 5)
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-15

for i in 0:5
    extrap1 = MidpointExtrapolation(pode, i)
    extrap2 = MidpointExtrapolationPODE(_get_v(pode), _get_f(pode), i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, q₀, qᵢ, p₀, pᵢ) == _midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, q₀, qᵢ, p₀, pᵢ) == _midpoint_extrapolation_pode!(_get_v(pode), _get_f(pode), t₀, t₁, q₀, qᵢ, p₀, pᵢ, i)
end

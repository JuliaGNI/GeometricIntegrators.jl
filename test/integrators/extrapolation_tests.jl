
using GeometricBase
using GeometricBase.Config
using GeometricEquations
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Tableaus
using GeometricProblems.HarmonicOscillator
using Test

using GeometricEquations: _get_v, _get_f, _get_v̄, _get_f̄
using GeometricIntegrators.Integrators: _euler_extrapolation_ode!, _hermite_extrapolation!,
                                        _midpoint_extrapolation_ode!, _midpoint_extrapolation_iode!, _midpoint_extrapolation_pode!

const Δt = 0.1
const t₋ = 0.0
const t₀ = 0.0
const t₁ = t₀ + Δt
const t₂ = t₁ + Δt

ode  = harmonic_oscillator_ode()
pode = harmonic_oscillator_pode()
iode = harmonic_oscillator_iode()


# Create ODE Solution Arrays

x₀ = ode.ics.q
x₁ = zero(x₀)
x₂ = zero(x₀)
xᵢ = zero(x₀)

ẋₚ = zero(x₀)
ẋ₀ = zero(x₀)
ẋ₁ = zero(x₀)
ẋ₂ = zero(x₀)
ẋₙ = zero(x₀)
ẋᵢ = zero(x₀)


# Compute Reference Solution for ODEs

ref_prev = integrate(similar(ode; tspan=(t₋,t₀), tstep=-Δt), TableauGauss(8), 1)
ref_next = integrate(similar(ode; tspan=(t₀,t₁), tstep=+Δt), TableauGauss(8), 1)

tₚ = ref_prev.t[end]
xₚ = ref_prev.q[end]

tₙ = ref_next.t[end]
xₙ = ref_next.q[end]

tᵢ = tₙ

functions(ode).v(tₚ, xₚ, ẋₚ)
functions(ode).v(t₀, x₀, ẋ₀)
functions(ode).v(tₙ, xₙ, ẋₙ)


# Hermite Extrapolation

evaluate!(HermiteExtrapolation(tₚ, t₀), xₚ, x₀, ẋₚ, ẋ₀, tᵢ, xᵢ, ẋᵢ)

# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, ẋₙ, ẋᵢ .- ẋₙ)

@test xᵢ ≈ xₙ atol=1E-5
@test ẋᵢ ≈ ẋₙ atol=1E-4

@test _hermite_extrapolation!(tₚ, t₀, xₚ, x₀, ẋₚ, ẋ₀, t₁, x₁) == xᵢ
@test _hermite_extrapolation!(tₚ, t₀, xₚ, x₀, ẋₚ, ẋ₀, t₁, x₁, ẋ₁) == (xᵢ, ẋᵢ)
@test _hermite_extrapolation!(functions(ode).v, tₚ, t₀, xₚ, x₀, t₁, x₁) == xᵢ
@test _hermite_extrapolation!(functions(ode).v, tₚ, t₀, xₚ, x₀, t₁, x₁, ẋ₁) == (xᵢ, ẋᵢ)


# Euler Extrapolation for ODEs

_euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 0)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-1

_euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 1)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-2

_euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 2)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-4

_euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 3)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-6

_euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 4)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-8

_euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 5)
# println(xᵢ, xₙ, xᵢ .- xₙ)
@test xᵢ ≈ xₙ atol=1E-10

for i in 0:5
    extrap1 = EulerExtrapolation(ode, i)
    extrap2 = EulerExtrapolationODE(functions(ode).v, i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, x₀, xᵢ) == _euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, x₀, xᵢ) == _euler_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, i)
end


# Midpoint Extrapolation for ODEs

_midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 0)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-4

_midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 1)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-8

_midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 2)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-12

_midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 3)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-15

_midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 4)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-16

_midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, 5)
# println(xᵢ, xₙ, qᵢ .- qₙ)
@test xᵢ ≈ xₙ atol=1E-15

for i in 0:5
    extrap1 = MidpointExtrapolation(ode, i)
    extrap2 = MidpointExtrapolationODE(functions(ode).v, i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, x₀, xᵢ) == _midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, x₀, xᵢ) == _midpoint_extrapolation_ode!(functions(ode).v, t₀, t₁, x₀, xᵢ, i)
end


# Create PODE Solution Arrays

q₀ = pode.ics.q
p₀ = pode.ics.p

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

ref_prev = integrate(similar(pode; tspan=(t₋,t₀), tstep=-Δt), PartitionedTableauGauss(8), 1)
ref_next = integrate(similar(pode; tspan=(t₀,t₁), tstep=+Δt), PartitionedTableauGauss(8), 1)

qₚ = ref_prev.q[end]
pₚ = ref_prev.p[end]

qₙ = ref_next.q[end]
pₙ = ref_next.p[end]

functions(pode).v(tₚ, qₚ, pₚ, q̇ₚ)
functions(pode).v(t₀, q₀, p₀, q̇₀)
functions(pode).v(tₙ, qₙ, pₙ, q̇ₙ)

functions(pode).f(tₚ, qₚ, pₚ, ṗₚ)
functions(pode).f(t₀, q₀, p₀, ṗ₀)
functions(pode).f(tₙ, qₙ, pₙ, ṗₙ)


# Midpoint Extrapolation for PODEs

_midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 0)
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-6
@test pᵢ ≈ pₙ atol=1E-4

_midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 1)
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-10
@test pᵢ ≈ pₙ atol=1E-8

_midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 2)
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-14
@test pᵢ ≈ pₙ atol=1E-12

_midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 3)
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

_midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 4)
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-16
@test pᵢ ≈ pₙ atol=1E-16

_midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 5)
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

for i in 0:5
    extrap1 = MidpointExtrapolation(pode, i)
    extrap2 = MidpointExtrapolationPODE(functions(pode).v, functions(pode).f, i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, q₀, qᵢ, p₀, pᵢ) == _midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, q₀, qᵢ, p₀, pᵢ) == _midpoint_extrapolation_pode!(functions(pode).v, functions(pode).f, t₀, t₁, q₀, qᵢ, p₀, pᵢ, i)
end


# Create IODE Solution Arrays

q₀ = iode.ics.q
p₀ = iode.ics.p

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


# Compute Reference Solution for IODEs

iode_prev = similar(iode; tspan=(t₋,t₀), tstep=-Δt)
iode_next = similar(iode; tspan=(t₀,t₁), tstep=+Δt)

ref_prev = integrate(iode_prev, IntegratorVPRKpSymmetric(iode_prev, TableauVPGLRK(8), tstep(iode_prev)), 1)
ref_next = integrate(iode_next, IntegratorVPRKpSymmetric(iode_next, TableauVPGLRK(8), tstep(iode_next)), 1)

qₚ = ref_prev.q[end]
pₚ = ref_prev.p[end]

qₙ = ref_next.q[end]
pₙ = ref_next.p[end]

functions(iode).v̄(tₚ, qₚ, q̇ₚ)
functions(iode).v̄(t₀, q₀, q̇₀)
functions(iode).v̄(tₙ, qₙ, q̇ₙ)

functions(iode).f̄(tₚ, qₚ, q̇ₚ, ṗₚ)
functions(iode).f̄(t₀, q₀, q̇₀, ṗ₀)
functions(iode).f̄(tₙ, qₙ, q̇ₙ, ṗₙ)



# Midpoint Extrapolation for IODEs

_midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 0)
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-4
@test pᵢ ≈ pₙ atol=1E-4

_midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 1)
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-8
@test pᵢ ≈ pₙ atol=1E-8

_midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 2)
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-12
@test pᵢ ≈ pₙ atol=1E-12

_midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 3)
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

_midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 4)
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-16
@test pᵢ ≈ pₙ atol=1E-16

_midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, 5)
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
@test qᵢ ≈ qₙ atol=1E-15
@test pᵢ ≈ pₙ atol=1E-16

for i in 0:5
    extrap1 = MidpointExtrapolation(iode, i)
    extrap2 = MidpointExtrapolationIODE(functions(iode).v̄, functions(iode).f̄, i)
    @test extrap1 == extrap2
    @test evaluate!(extrap1, t₀, t₁, q₀, qᵢ, p₀, pᵢ) == _midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, i)
    @test evaluate!(extrap2, t₀, t₁, q₀, qᵢ, p₀, pᵢ) == _midpoint_extrapolation_iode!(functions(iode).v̄, functions(iode).f̄, t₀, t₁, q₀, qᵢ, p₀, pᵢ, i)
end

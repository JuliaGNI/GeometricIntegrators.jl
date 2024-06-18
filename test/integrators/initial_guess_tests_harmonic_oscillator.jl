using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: Δt
using Test

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()

dae  = daeproblem()
pdae = pdaeproblem()
hdae = hdaeproblem()
idae = idaeproblem()
ldae = ldaeproblem()


# ODE Reference Solution

ode_prev = similar(ode; tspan=(tspan(ode)[begin], tspan(ode)[begin]-tstep(ode)), tstep=-tstep(ode))
ode_next = similar(ode; tspan=(tspan(ode)[begin], tspan(ode)[begin]+tstep(ode)), tstep=+tstep(ode))

t₀ = initial_conditions(ode).t
q₀ = initial_conditions(ode).q
v₀ = zero(q₀)

tₚ = tspan(ode_prev)[end]
qₚ = zero(q₀)
vₚ = zero(v₀)

tₙ = tspan(ode_next)[end]
qₙ = zero(q₀)
vₙ = zero(v₀)

hist = (t = [t₀], q = [q₀])
prev = (t = tₚ, q = qₚ)
next = (t = tₙ, q = qₙ)
extrapolate!(prev, hist, ode_prev, MidpointExtrapolation(5))
extrapolate!(next, hist, ode_next, MidpointExtrapolation(5))

equation(ode).v(vₚ, tₚ, qₚ, parameters(ode))
equation(ode).v(vₙ, tₙ, qₙ, parameters(ode))


# ODE Initial Guess

t₀ = initial_conditions(ode).t
q₀ = initial_conditions(ode).q
v₀ = zero(q₀)

t₁ = t₀ - Δt
q₁ = zero(q₀)
v₁ = zero(v₀)

t₂ = t₀ + Δt
q₂ = zero(q₀)
v₂ = zero(v₀)

equation(ode).v(v₀, t₀, q₀, parameters(ode))

initialguess!(t₀, q₀, t₁, q₁, v₁, ode, MidpointExtrapolation(4))
initialguess!(t₁, q₁, v₁, t₀, q₀, v₀, t₂, q₂, v₂, HermiteExtrapolation())

# println("ODE Initial Guess")
# println("Δq = $(q₁ .- qₚ)")
# println("Δv = $(v₁ .- vₚ)")
# println("Δq = $(q₂ .- qₙ)")
# println("Δv = $(v₂ .- vₙ)")
# println()

@test q₁ ≈ qₚ atol=1E-15
@test v₁ ≈ vₚ atol=1E-15

@test q₂ ≈ qₙ atol=1E-5
@test v₂ ≈ vₙ atol=1E-4


# PODE Reference Solution

pode_prev = similar(pode; tspan=(tspan(pode)[begin], tspan(pode)[begin]-tstep(pode)), tstep=-tstep(pode))
pode_next = similar(pode; tspan=(tspan(pode)[begin], tspan(pode)[begin]+tstep(pode)), tstep=+tstep(pode))

t₀ = initial_conditions(pode).t
q₀ = initial_conditions(pode).q
p₀ = initial_conditions(pode).p
v₀ = zero(q₀)
f₀ = zero(p₀)

tₚ = tspan(pode_prev)[end]
qₚ = zero(q₀)
pₚ = zero(p₀)
vₚ = zero(v₀)
fₚ = zero(p₀)

tₙ = tspan(pode_next)[end]
qₙ = zero(q₀)
pₙ = zero(p₀)
vₙ = zero(v₀)
fₙ = zero(f₀)

hist = (t = [t₀], q = [q₀], p = [p₀])
prev = (t = tₚ, q = qₚ, p = pₚ)
next = (t = tₙ, q = qₙ, p = pₙ)
extrapolate!(prev, hist, pode_prev, MidpointExtrapolation(5))
extrapolate!(next, hist, pode_next, MidpointExtrapolation(5))

equation(pode).v(vₚ, tₚ, qₚ, pₚ, parameters(pode))
equation(pode).v(vₙ, tₙ, qₙ, pₙ, parameters(pode))

equation(pode).f(fₚ, tₚ, qₚ, pₚ, parameters(pode))
equation(pode).f(fₙ, tₙ, qₙ, pₙ, parameters(pode))


# PODE Initial Guess

t₀ = initial_conditions(pode).t
q₀ = initial_conditions(pode).q
p₀ = initial_conditions(pode).p
v₀ = zero(q₀)
f₀ = zero(p₀)

t₁ = t₀ - Δt
q₁ = zero(q₀)
p₁ = zero(p₀)
v₁ = zero(v₀)
f₁ = zero(f₀)

t₂ = t₀ + Δt
q₂ = zero(q₀)
p₂ = zero(p₀)
v₂ = zero(v₀)
f₂ = zero(f₀)

equation(pode).v(v₀, t₀, q₀, p₀, parameters(pode))
equation(pode).f(f₀, t₀, q₀, p₀, parameters(pode))

initialguess!(t₀, q₀, p₀, t₁, q₁, p₁, v₁, f₁, pode, MidpointExtrapolation(4))
initialguess!(t₁, q₁, p₁, v₁, f₁, t₀, q₀, p₀, v₀, f₀, t₂, q₂, p₂, v₂, f₂, HermiteExtrapolation())

# println("PODE Initial Guess")
# println("Δq = $(q₁ .- qₚ)")
# println("Δp = $(p₁ .- pₚ)")
# println("Δv = $(v₁ .- vₚ)")
# println("Δf = $(f₁ .- fₚ)")
# println("Δq = $(q₂ .- qₙ)")
# println("Δp = $(p₂ .- pₙ)")
# println("Δv = $(v₂ .- vₙ)")
# println("Δf = $(f₂ .- fₙ)")
# println()

@test q₁ ≈ qₚ atol=1E-15
@test p₁ ≈ pₚ atol=1E-15
@test v₁ ≈ vₚ atol=1E-15
@test f₁ ≈ fₚ atol=1E-15

@test q₂ ≈ qₙ atol=1E-5
@test p₂ ≈ pₙ atol=1E-7
@test v₂ ≈ vₙ atol=1E-4
@test f₂ ≈ fₙ atol=1E-6


# IODE Initial Guess

t₀ = initial_conditions(iode).t
q₀ = initial_conditions(iode).q
p₀ = initial_conditions(iode).p
v₀ = zero(q₀)
f₀ = zero(p₀)

t₁ = t₀ - Δt
q₁ = zero(q₀)
p₁ = zero(p₀)
v₁ = zero(v₀)
f₁ = zero(f₀)

t₂ = t₀ + Δt
q₂ = zero(q₀)
p₂ = zero(p₀)
v₂ = zero(v₀)
f₂ = zero(f₀)

equation(iode).v̄(v₀, t₀, q₀, p₀, parameters(iode))
equation(iode).f̄(f₀, t₀, q₀, v₀, parameters(iode))

initialguess!(t₀, q₀, p₀, t₁, q₁, p₁, v₁, f₁, iode, MidpointExtrapolation(4))
initialguess!(t₁, q₁, p₁, v₁, f₁, t₀, q₀, p₀, v₀, f₀, t₂, q₂, p₂, v₂, f₂, HermiteExtrapolation())

# println("IODE Initial Guess")
# println("Δq = $(q₁ .- qₚ)")
# println("Δp = $(p₁ .- pₚ)")
# println("Δv = $(v₁ .- vₚ)")
# println("Δf = $(f₁ .- fₚ)")
# println("Δq = $(q₂ .- qₙ)")
# println("Δp = $(p₂ .- pₙ)")
# println("Δv = $(v₂ .- vₙ)")
# println("Δf = $(f₂ .- fₙ)")
# println()

@test q₁ ≈ qₚ atol=1E-15
@test p₁ ≈ pₚ atol=1E-15
@test v₁ ≈ vₚ atol=1E-15
@test f₁ ≈ fₚ atol=1E-15

@test q₂ ≈ qₙ atol=1E-5
@test p₂ ≈ pₙ atol=1E-7
@test v₂ ≈ vₙ atol=1E-4
@test f₂ ≈ fₙ atol=1E-6

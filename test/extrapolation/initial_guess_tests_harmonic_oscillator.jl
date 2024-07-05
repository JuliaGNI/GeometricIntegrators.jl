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
equation(ode).v(v₀, t₀, q₀, parameters(ode))

tₚ = tspan(ode_prev)[end]
qₚ = zero(q₀)
vₚ = zero(v₀)

tₙ = tspan(ode_next)[end]
qₙ = zero(q₀)
vₙ = zero(v₀)

hist = (t = [t₀], q = [q₀], v = [v₀])
prev = (t = tₚ, q = qₚ, v = vₚ)
next = (t = tₙ, q = qₙ, v = vₙ)
solutionstep!(prev, hist, ode_prev, MidpointExtrapolation(5))
solutionstep!(next, hist, ode_next, MidpointExtrapolation(5))


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

hist = (t = [t₀], q = [q₀], v = [v₀])
sol1 = (t = t₁, q = q₁, v = v₁)
solutionstep!(sol1, hist, ode, MidpointExtrapolation(4))

hist = (t = [t₁, t₀], q = [q₁, q₀], v = [v₁, v₀])
sol2 = (t = t₂, q = q₂, v = v₂)
solutionstep!(sol2, hist, ode, HermiteExtrapolation())

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

hist = (t = [t₀], q = [q₀], p = [p₀], v = [v₀], f = [f₀])
prev = (t = tₚ, q = qₚ, p = pₚ, v = vₚ, f = fₚ)
next = (t = tₙ, q = qₙ, p = pₙ, v = vₙ, f = fₙ)
solutionstep!(prev, hist, pode_prev, MidpointExtrapolation(5))
solutionstep!(next, hist, pode_next, MidpointExtrapolation(5))


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

hist = (t = [t₀], q = [q₀], p = [p₀], v = [v₀], f = [f₀])
sol1 = (t = t₁, q = q₁, p = p₁, v = v₁, f = f₁)
solutionstep!(sol1, hist, pode, MidpointExtrapolation(4))

hist = (t = [t₁, t₀], q = [q₁, q₀], p = [p₁, p₀], v = [v₁, v₀], f = [f₁, f₀])
sol2 = (t = t₂, q = q₂, p = p₂, v = v₂, f = f₂)
solutionstep!(sol2, hist, pode, HermiteExtrapolation())

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

initialguess(iode).v(v₀, t₀, q₀, p₀, parameters(iode))
initialguess(iode).f(f₀, t₀, q₀, v₀, parameters(iode))

hist = (t = [t₀], q = [q₀], p = [p₀], v = [v₀], f = [f₀])
sol1 = (t = t₁, q = q₁, p = p₁, v = v₁, f = f₁)
solutionstep!(sol1, hist, iode, MidpointExtrapolation(4))

hist = (t = [t₁, t₀], q = [q₁, q₀], p = [p₁, p₀], v = [v₁, v₀], f = [f₁, f₀])
sol2 = (t = t₂, q = q₂, p = p₂, v = v₂, f = f₂)
solutionstep!(sol2, hist, iode, HermiteExtrapolation())

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

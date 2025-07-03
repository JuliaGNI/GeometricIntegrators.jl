import GeometricEquations: DELEProblem

abstract type DELEQuadrature end

struct Midpoint <: DELEQuadrature end
struct Trapezoidal <: DELEQuadrature end


function dele_initial_data(prob::LODEProblem, predictor::LODEMethod)
    t₀ = timespan(prob)[begin]
    t₁ = t₀ - timestep(prob)
    sol = integrate(similar(prob, (t₀, t₁), -timestep(prob)), predictor)
    return (sol[1].q, sol[0].q)
end


function dele_midpoint_Ld(t₀, t₁, q₀, q₁, params, prob)
    h = (t₁ - t₀)
    t = (t₀ + t₁) / 2
    q = (q₀ .+ q₁) ./ 2
    v = (q₁ .- q₀) ./ h

    return h * functions(prob).lagrangian(t, q, v, params)
end

function dele_midpoint_D1Ld(d, t₀, t₁, q₀, q₁, params, prob)
    h = (t₁ - t₀)
    t = (t₀ + t₁) / 2
    q = (q₀ .+ q₁) ./ 2
    v = (q₁ .- q₀) ./ h
    p = zero(q)
    f = zero(q)

    functions(prob).ϑ(p, t, q, v, params)
    functions(prob).f(f, t, q, v, params)

    d .= h .* f ./ 2 .- p

    return nothing
end

function dele_midpoint_D2Ld(d, t₀, t₁, q₀, q₁, params, prob)
    h = (t₁ - t₀)
    t = (t₀ + t₁) / 2
    q = (q₀ .+ q₁) ./ 2
    v = (q₁ .- q₀) ./ h
    p = zero(q)
    f = zero(q)

    functions(prob).ϑ(p, t, q, v, params)
    functions(prob).f(f, t, q, v, params)

    d .= h .* f ./ 2 .+ p

    return nothing
end

function GeometricEquations.DELEProblem(prob::LODEProblem, ::Midpoint, predictor=VPRKGauss(8))
    q₀, q₁ = dele_initial_data(prob, predictor)

    DELEProblem(
        (t₀, t₁, q₀, q₁, params) -> dele_midpoint_Ld(t₀, t₁, q₀, q₁, params, prob),
        (d, t₀, t₁, q₀, q₁, params) -> dele_midpoint_D1Ld(d, t₀, t₁, q₀, q₁, params, prob),
        (d, t₀, t₁, q₀, q₁, params) -> dele_midpoint_D2Ld(d, t₀, t₁, q₀, q₁, params, prob),
        timespan(prob), timestep(prob), q₀, q₁; invariants=invariants(prob), parameters=parameters(prob))
end



function dele_trapezoidal_Ld(t₀, t₁, q₀, q₁, params, prob)
    h = (t₁ - t₀)
    v = (q₁ .- q₀) ./ h

    return h * (functions(prob).lagrangian(t₀, q₀, v, params) + functions(prob).lagrangian(t₁, q₁, v, params)) / 2
end

function dele_trapezoidal_D1Ld(d, t₀, t₁, q₀, q₁, params, prob)
    h = (t₁ - t₀)
    v = (q₁ .- q₀) ./ h

    p₀ = zero(v)
    p₁ = zero(v)
    f₀ = zero(v)

    functions(prob).ϑ(p₀, t₀, q₀, v, params)
    functions(prob).ϑ(p₁, t₁, q₁, v, params)
    functions(prob).f(f₀, t₀, q₀, v, params)

    d .= h .* f₀ ./ 2 .- p₀ ./ 2 .- p₁ ./ 2

    return nothing
end

function dele_trapezoidal_D2Ld(d, t₀, t₁, q₀, q₁, params, prob)
    h = (t₁ - t₀)
    v = (q₁ .- q₀) ./ h

    p₀ = zero(v)
    p₁ = zero(v)
    f₁ = zero(v)

    functions(prob).ϑ(p₀, t₀, q₀, v, params)
    functions(prob).ϑ(p₁, t₁, q₁, v, params)
    functions(prob).f(f₁, t₁, q₁, v, params)

    d .= h .* f₁ ./ 2 .+ p₀ ./ 2 .+ p₁ ./ 2

    return nothing
end

function GeometricEquations.DELEProblem(prob::LODEProblem, ::Trapezoidal, predictor=VPRKGauss(8))
    q₀, q₁ = dele_initial_data(prob, predictor)

    DELEProblem(
        (t₀, t₁, q₀, q₁, params) -> dele_trapezoidal_Ld(t₀, t₁, q₀, q₁, params, prob),
        (d, t₀, t₁, q₀, q₁, params) -> dele_trapezoidal_D1Ld(d, t₀, t₁, q₀, q₁, params, prob),
        (d, t₀, t₁, q₀, q₁, params) -> dele_trapezoidal_D2Ld(d, t₀, t₁, q₀, q₁, params, prob),
        timespan(prob), timestep(prob), q₀, q₁; invariants=invariants(prob), parameters=parameters(prob))
end

@doc raw"""
Projected Gauss-Legendre Runge-Kutta method: an energy-preserving Gauss collocation
method for [`ODE`](@ref)s carrying an invariant `h`.

The method uses the one-parameter family of tableaus of [`CoefficientsPGLRK`](@ref),
```math
a(\lambda) = a + \lambda A ,
\qquad
A = P W Q ,
```
where ``a`` is the ``s``-stage Gauss tableau in W-transformed form and ``W`` is skew.
The stages
```math
\begin{aligned}
Y_{n,i} &= \sum \limits_{j=1}^{s} a_{ij} (\lambda) \, V_{n,j} , &
Q_{n,i} &= q_{n} + h \, Y_{n,i} , &
V_{n,i} &= v (t_{n,i}, Q_{n,i}) , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} ,
\end{aligned}
```
are solved for each trial ``\lambda``, and ``\lambda`` is then determined by a
bisection on the scalar energy-conservation condition
```math
H (q_{n+1} (\lambda)) = H (q_{0}) .
```

Because ``B A`` is skew for ``B = \mathrm{diag}(b)``, the perturbed *tableau* satisfies
the symplecticity condition ``b_{i} a_{ij} + b_{j} a_{ji} = b_{i} b_{j}`` for every
*fixed* ``\lambda``, and because ``b^{T} A = 0`` and ``A \mathbb{1} = 0`` it also retains
the quadrature and row-sum conditions, hence the order ``2s``. Both identities are
asserted in `test/methods/pglrk_coefficients_tests.jl`.

That does **not** make the method symplectic. ``\lambda`` is determined from ``q_{n}`` at
every step, so the step map is a member of the ``\lambda``-family composed with a
state-dependent choice of parameter, and its Jacobian carries an additional rank-one term
``(\partial \Psi / \partial \lambda) (\partial \lambda / \partial q)^{T}``. Exact energy
conservation and symplecticity are in any case mutually exclusive for a general
Hamiltonian system unless the method reproduces the exact flow (Ge & Marsden, 1988).
Measured on a nonlinear pendulum at ``\Delta t = 0.8``, the defect
``\vert J^{T} \Omega J - \Omega \vert`` tends to ``3.1 \times 10^{-7}`` as the
finite-difference step is refined, against ``7.3 \times 10^{-12}`` for `Gauss(3)` — the
latter being finite-difference noise rather than a defect. Accordingly `issymplectic`
returns `missing`.

The reference energy is that of the **initial condition of the problem**, not of the
previous step, so the absolute energy error stays bounded instead of performing a random
walk.

Note that ``\lambda`` is sought in the bracket ``[-\lambda_{\max}, +\lambda_{\max}]``
with ``\lambda_{\max} = 2 \cdot 10^{-s}`` by default. If the energy residual does not
change sign across that bracket there is no root, and the method falls back to
``\lambda = 0``, i.e. to the plain Gauss method, rather than applying a large
perturbation; the number of steps on which that happens is counted in the cache's
`nfallback`.

!!! note
    The scalar ``\lambda`` solve uses `SimpleSolvers.bisection`, which SimpleSolvers does
    not export, so a compat bump there may break this method.

### Reference

Luigi Brugnano, Felice Iavernaro and Donato Trigiante.
Energy- and quadratic invariants-preserving integrators based upon Gauss collocation
formulae.
SIAM Journal on Numerical Analysis 50(6), pp. 2897-2916, 2012.
doi: [10.1137/110856617](https://doi.org/10.1137/110856617).

### Constructors

```
PGLRK(coefficients::CoefficientsPGLRK; λmax = 2 / 10^s)
PGLRK(s::Int, T = Float64; λmax = 2 / 10^s)
```
"""
struct PGLRK{CT<:CoefficientsPGLRK,T} <: ODEMethod
    coefficients::CT
    λmax::T

    # `T(10)^s` rather than `10^s`: the latter overflows `Int64` for s ≥ 19
    function PGLRK(coefficients::CoefficientsPGLRK{T}; λmax=T(2) / T(10)^RungeKutta.nstages(coefficients)) where {T}
        @assert λmax > 0 "λmax must be positive"
        new{typeof(coefficients),T}(coefficients, λmax)
    end
end

PGLRK(s::Int, ::Type{T}=Float64; kwargs...) where {T} = PGLRK(CoefficientsPGLRK(T, s); kwargs...)

GeometricBase.tableau(method::PGLRK, args...; kwargs...) = method.coefficients
GeometricBase.order(method::PGLRK) = RungeKutta.order(tableau(method))
GeometricBase.order(::Type{<:PGLRK}) = "2s"

@inline nstages(method::PGLRK) = RungeKutta.nstages(tableau(method))
@inline eachstage(method::PGLRK) = RungeKutta.eachstage(tableau(method))

isexplicit(::Union{PGLRK,Type{<:PGLRK}}) = false
isimplicit(::Union{PGLRK,Type{<:PGLRK}}) = true
issymmetric(::Union{PGLRK,Type{<:PGLRK}}) = true

# The *tableau* a(λ) is symplectic for every fixed λ (`B A` is skew), but λ is solved per
# step from the energy condition, so the step map is that tableau composed with a
# state-dependent parameter choice and its Jacobian carries an extra rank-one term.
# Ge & Marsden (1988) also rule out a method that conserves H exactly *and* is symplectic
# unless it reproduces the exact flow. Measured on a nonlinear pendulum at Δt = 0.8:
# |JᵀΩJ − Ω| → 3.1E-7 as the finite-difference step is refined, against 7.3E-12 for
# Gauss(3) — i.e. Gauss's residual is FD noise and PGLRK's is a genuine defect.
issymplectic(::Union{PGLRK,Type{<:PGLRK}}) = missing
isenergypreserving(::Union{PGLRK,Type{<:PGLRK}}) = true

default_solver(::PGLRK) = Newton()
default_iguess(::PGLRK) = HermiteExtrapolation()

solversize(method::PGLRK, problem::AbstractProblemODE) =
    length(vec(initial_conditions(problem).q)) * nstages(method)


function Base.show(io::IO, int::GeometricIntegrator{<:PGLRK})
    print(io, "\nProjected Gauss-Legendre Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int)))")
end


@doc raw"""
Projected Gauss-Legendre Runge-Kutta integrator cache.

Mutable, because the projection parameter `λ`, the reference energy `h₀` and the
working tableau `ā = a + λ A` are updated in the course of a step.

### Fields

* `x`: nonlinear solver solution vector, holding `Y`
* `x₀`: λ-independent starting point for the stage solves
* `q`: candidate solution ``q_{n+1}`` for the current `λ`
* `q̃`, `ṽ`: temporaries for the initial guess
* `Q`, `V`, `Y`: internal stages
* `ā`: working tableau ``a + \lambda A``
* `λ`: current value of the projection parameter
* `h`, `h₀`: current and reference value of the invariant
* `nfallback`: number of steps on which no root was found in ``[-\lambda_{\max}, +\lambda_{\max}]``
  and the method fell back to ``\lambda = 0``, i.e. to the plain Gauss method
"""
mutable struct PGLRKCache{ST,D,S} <: ODEIntegratorCache{ST}
    x::Vector{ST}
    x₀::Vector{ST}

    q::Vector{ST}
    q̃::Vector{ST}
    ṽ::Vector{ST}

    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    Y::Vector{Vector{ST}}

    ā::Matrix{ST}

    λ::ST
    h::ST
    h₀::ST
    nfallback::Int

    function PGLRKCache{ST,D,S}(coefficients) where {ST,D,S}
        x = zeros(ST, D * S)
        x₀ = zeros(ST, D * S)

        q = zeros(ST, D)
        q̃ = zeros(ST, D)
        ṽ = zeros(ST, D)

        Q = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)

        # Start from the unperturbed Gauss tableau. This MUST be a fresh array:
        # `convert(Matrix{ST}, x)` returns `x` itself when it is already a `Matrix{ST}`,
        # which would alias `ā` to the method's own tableau, so that every
        # `update_tableau!` would accumulate into — and permanently corrupt — it.
        ā = Matrix{ST}(coefficients.a)

        new(x, x₀, q, q̃, ṽ, Q, V, Y, ā, zero(ST), zero(ST), zero(ST), 0)
    end
end

nlsolution(cache::PGLRKCache) = cache.x

function Cache{ST}(problem::AbstractProblemODE, method::PGLRK; kwargs...) where {ST}
    @assert hasinvariants(problem) "PGLRK requires a problem with an invariant `h`"
    @assert :h ∈ keys(invariants(problem)) "PGLRK requires an invariant named `h`"
    D = length(vec(initial_conditions(problem).q))
    PGLRKCache{ST,D,nstages(method)}(tableau(method); kwargs...)
end

@inline function CacheType(ST, problem::AbstractProblemODE, method::PGLRK)
    D = length(vec(initial_conditions(problem).q))
    PGLRKCache{ST,D,nstages(method)}
end


function internal_variables(method::PGLRK, problem::AbstractProblemODE{DT,TT}) where {DT,TT}
    S = nstages(method)
    D = length(vec(initial_conditions(problem).q))

    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    Y = create_internal_stage_vector(DT, D, S)

    (Q=Q, V=V, Y=Y)
end

function copy_internal_variables!(solstep::SolutionStep, cache::PGLRKCache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :Y) && copyto!(internal(solstep).Y, cache.Y)
end


"""
Set the projection parameter `λ` and recompute the working tableau `ā = a + λ A`.
Both live in the `DT` cache, since `λ` is real-valued and constant with respect to
the nonlinear solver unknowns.
"""
function update_tableau!(int::GeometricIntegrator{<:PGLRK}, λ)
    cache(int).λ = λ
    getTableauPGLRK(tableau(int), λ, cache(int).ā)
    return cache(int).ā
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:PGLRK})
    local x = nlsolution(int)
    local D = length(cache(int).q̃)

    # compute initial guess for the internal stages
    for i in eachstage(int)
        soltmp = (
            t=history[1].t + timestep(int) * tableau(int).c[i],
            q=cache(int).Q[i],
            q̇=cache(int).V[i],
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))
    end

    # assemble the solution vector, which holds Y = a V (no Δt factor)
    for i in eachstage(int)
        offset = D * (i - 1)
        for k in 1:D
            x[offset+k] = 0
            for j in eachstage(int)
                x[offset+k] += tableau(int).a[i, j] * cache(int).V[j][k]
            end
        end
    end

    # remember it: every trial λ must start the stage solve from the same point, or
    # `energy_residual!` becomes path-dependent and the bisection wanders
    copyto!(cache(int).x₀, x)
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:PGLRK}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)

    # copy x to Y
    for i in eachindex(C.Y)
        for k in eachindex(C.Y[i])
            C.Y[i][k] = x[D*(i-1)+k]
        end
    end

    # compute Q = q + Δt Y
    for i in eachindex(C.Q)
        for k in eachindex(C.Q[i])
            C.Q[i][k] = sol.q[k] + timestep(int) * C.Y[i][k]
        end
    end

    # compute V = v(Q)
    for i in eachindex(C.Q, C.V)
        tᵢ = sol.t + timestep(int) * (tableau(int).c[i] - 1)
        equations(int).v(C.V[i], tᵢ, C.Q[i], params)
    end

    # compute the candidate solution q = q + Δt b V
    for k in eachindex(C.q)
        y = zero(ST)
        for j in eachindex(C.V)
            y += tableau(int).b[j] * C.V[j][k]
        end
        C.q[k] = sol.q[k] + timestep(int) * y
    end
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:PGLRK}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)

    # the λ-dependent tableau is real-valued and is read from the DT cache
    local ā = cache(int).ā

    # compute b = Y - a(λ) V. Note that `Y` here carries *no* Δt factor (unlike `FLRK`,
    # where the unknown is Y = Δt a V); the factor sits in `components!` instead.
    for i in eachindex(C.Y)
        for k in eachindex(C.Y[i])
            y = zero(ST)
            for j in eachindex(C.V)
                y += ā[i, j] * C.V[j][k]
            end
            b[D*(i-1)+k] = C.Y[i][k] - y
        end
    end
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:PGLRK}) where {ST}
    @assert axes(x) == axes(b)
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end


"""
Solve the stage equations for a trial `λ` and return the energy residual
`h₀ - h(q_{n+1}(λ))`.

The stage solve is restarted from the λ-independent initial guess `x₀` on every call.
Warm-starting from the previous trial `λ` instead makes the returned residual depend
on the sequence of trials rather than on `λ` alone, which breaks the bisection: the
bracket then collapses to a denormal instead of converging on the root.
"""
function energy_residual!(λ, (sol, params, int))
    update_tableau!(int, λ)
    copyto!(nlsolution(int), cache(int).x₀)
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))
    components!(nlsolution(int), sol, params, int)
    cache(int).h = invariants(problem(int)).h(sol.t, cache(int).q, params)
    return cache(int).h₀ - cache(int).h
end


"""
Determine the projection parameter `λ` by bisection on the energy residual.

Falls back to `λ = 0` (the plain Gauss method) when the residual already vanishes, or
when it does not change sign across the bracket — in the latter case
`SimpleSolvers.bisection` returns an endpoint, i.e. the largest admissible perturbation
rather than a root. Such steps are counted in `cache(int).nfallback`.
"""
function solve_λ!(sol, params, int::GeometricIntegrator{<:PGLRK,<:AbstractProblemODE}, DT)
    local C = cache(int)
    local λmax = DT(method(int).λmax)
    local ftol = max(abs(C.h₀), one(DT)) * 8eps(DT)
    local args = (sol, params, int)

    # the unperturbed method may already conserve the energy to tolerance
    abs(energy_residual!(zero(DT), args)) ≤ ftol && return zero(DT)

    # Do *not* pre-probe ±λmax here: every probe is a full nonlinear stage solve, and
    # `bisection` evaluates both endpoints itself. When the bracket shows no sign change
    # it returns the endpoint with the smaller |f|, which is the largest admissible
    # perturbation rather than a root — detected and rejected below.
    λ = SimpleSolvers.bisection(energy_residual!, -λmax, +λmax, args,
        Options(DT; f_abstol=ftol, x_suctol=8eps(DT), verbosity=0))

    if abs(λ) ≥ λmax
        @debug "PGLRK: no energy-conserving λ in [±$(λmax)]; falling back to plain Gauss."
        C.nfallback += 1
        λ = zero(DT)
    end

    # leave the cache consistent with the accepted λ
    energy_residual!(λ, args)

    return λ
end


function update!(sol, params, int::GeometricIntegrator{<:PGLRK}, DT)
    # the accepted candidate is exactly the q whose energy was matched
    sol.q .= cache(int, DT).q
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:PGLRK,<:AbstractProblemODE})
    local DT = eltype(cache(int).q)

    # Reference energy, read from the problem's initial condition rather than latched on
    # the first step. Recomputing it costs one scalar evaluation against the ~10 nonlinear
    # stage solves below, and it removes a piece of mutable state whose invariant ("h₀
    # belongs to the initial condition this integrator was built for") would otherwise
    # hold only by construction. Evaluated with the step's own `params` so that a caller
    # varying parameters between runs stays correct. Note `ic.t` is a zero-dimensional
    # array and has to be indexed.
    local ic = initial_conditions(problem(int))
    cache(int).h₀ = invariants(problem(int)).h(ic.t[], ic.q, params)

    # determine λ from the energy-conservation condition; this performs the stage
    # solves and leaves the cache consistent with the accepted λ
    solve_λ!(sol, params, int, DT)

    # compute the final update
    update!(sol, params, int, DT)
end

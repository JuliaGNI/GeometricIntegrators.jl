@doc raw"""
Abstract supertype of the Discontinuous Galerkin Variational Integrators.

DGVIs discretise the action of a **fully degenerate** Lagrangian
```math
L(q, \dot{q}) = \vartheta (q) \cdot \dot{q} - H(q)
```
by a piecewise-polynomial trajectory that is *not* required to be continuous across
the interval boundaries, plus a numerical flux that integrates the Lagrangian across
the resulting jumps. Within each interval,
```math
q_h(t) \vert_{(t_{n}, t_{n+1})} = \sum \limits_{i=1}^{S} x_{n,i} \, \varphi_{i} \bigg( \frac{t - t_{n}}{t_{n+1} - t_{n}} \bigg) ,
```
and with a quadrature rule of ``R`` nodes ``c_i`` and weights ``b_i`` the mass and
derivative matrices and the boundary reconstruction vectors are
```math
m_{ij} = \varphi_j (c_i) ,
\qquad
a_{ij} = \varphi_j' (c_i) ,
\qquad
r^{+}_{j} = \varphi_j (0) ,
\qquad
r^{-}_{j} = \varphi_j (1) ,
```
so that
```math
Q_{n,i} = m_{ij} x_{n,j} ,
\qquad
V_{n,i} = \frac{a_{ij}}{h} x_{n,j} ,
\qquad
q_{n}^{+} = r^{+}_{j} x_{n,j} ,
\qquad
q_{n+1}^{-} = r^{-}_{j} x_{n,j} .
```
All variants share the interior contribution to the discrete Euler-Lagrange equations,
```math
0 = \sum \limits_{i=1}^{R} b_i \big[ h \, m_{ij} \, F_{n,i} + a_{ij} \, P_{n,i} \big] + (\text{flux terms}) ,
```
with ``P_{n,i} = \vartheta (Q_{n,i})`` and ``F_{n,i} = f (Q_{n,i}, V_{n,i})``, and differ
only in the flux and in the equation that closes the step. See the individual method
docstrings and `docs/src/integrators/dgvi.md`.

The five variants are [`DGVI`](@ref), [`DGVIPI`](@ref), [`DGVIP0`](@ref),
[`DGVIP1`](@ref) and [`DGVIEXP`](@ref).

!!! warning
    These methods require a degenerate Lagrangian: `ϑ` is evaluated as `ϑ(θ, t, q, q)`,
    which is only meaningful when it does not depend on the velocity. On a regular
    Lagrangian the closure equation degenerates — for `DGVI` on the harmonic oscillator
    it collapses to `q - p`, which is independent of every unknown and yields a singular
    Jacobian. Use `GeometricProblems.LotkaVolterra2d.iodeproblem_dg`,
    `GeometricProblems.HarmonicOscillator.degenerate_lodeproblem` or similar.
"""
abstract type DGVIMethod{T, S, R, F} <: LODEMethod end

isexplicit(::Union{DGVIMethod, Type{<:DGVIMethod}}) = false
isimplicit(::Union{DGVIMethod, Type{<:DGVIMethod}}) = true
issymmetric(::Union{DGVIMethod, Type{<:DGVIMethod}}) = missing
issymplectic(::Union{DGVIMethod, Type{<:DGVIMethod}}) = missing
isiodemethod(::Union{DGVIMethod, Type{<:DGVIMethod}}) = true

default_solver(::DGVIMethod) = Newton()
default_iguess(::DGVIMethod) = HermiteExtrapolation()

basis(method::DGVIMethod) = method.basis
quadrature(method::DGVIMethod) = method.quadrature

nbasis(::DGVIMethod{T, S, R, F}) where {T, S, R, F} = S
nnodes(::DGVIMethod{T, S, R, F}) where {T, S, R, F} = R

"""
Number of quadrature nodes of the jump path integral; zero for every variant except
[`DGVIPI`](@ref).
"""
njump(::DGVIMethod{T, S, R, F}) where {T, S, R, F} = F

"""
Number of degrees of freedom the step adds on top of the `S` basis coefficients.
"""
nclosure(::DGVIMethod) = 1

function solversize(method::DGVIMethod, problem::AbstractProblemIODE)
    length(vec(initial_conditions(problem).q)) * (nbasis(method) + nclosure(method))
end

"""
Compute the shared DGVI coefficient block `(b, c, x, m, a, r⁻, r⁺)` from a basis and a
quadrature rule.
"""
function dgvi_coefficients(basis::Basis{T}, quadrature::QuadratureRule{T}) where {T}
    R = QuadratureRules.nnodes(quadrature)
    S = CompactBasisFunctions.nbasis(basis)

    quad_nodes = QuadratureRules.nodes(quadrature)
    quad_weights = QuadratureRules.weights(quadrature)

    m = zeros(T, R, S)
    a = zeros(T, R, S)
    r⁻ = zeros(T, S)
    r⁺ = zeros(T, S)

    for i in eachindex(basis)
        r⁺[i] = basis[zero(T), i]
        r⁻[i] = basis[one(T), i]
        for j in eachindex(quad_nodes)
            m[j, i] = basis[quad_nodes[j], i]
            a[j, i] = basis'[quad_nodes[j], i]
        end
    end

    return (
        b = SVector{R, T}(quad_weights),
        c = SVector{R, T}(quad_nodes),
        x = SVector{S, T}(CompactBasisFunctions.grid(basis)),
        m = SMatrix{R, S, T, R * S}(m),
        a = SMatrix{R, S, T, R * S}(a),
        r⁻ = SVector{S, T}(r⁻),
        r⁺ = SVector{S, T}(r⁺)
    )
end

"""
Compute the jump path-integral coefficient block
`(β, γ, μ⁻, μ⁺, α⁻, α⁺, ρ⁻, ρ⁺)` of a [`Discontinuity`](@ref).

`μ∓` and `α∓` are the "mass" and "derivative" vectors of the connecting path,
`μ^∓_i = φ^∓(γ_i)` and `α^∓_i = dφ^∓/dτ(γ_i)`, and `ρ∓ = φ^∓(1/2)` reconstruct the
nodal value `q_n = Φ(1/2; q_n^-, q_n^+)`.
"""
function dgvi_jump_coefficients(jump::Discontinuity{T}) where {T}
    F = QuadratureRules.nnodes(jump.quadrature)

    γ = QuadratureRules.nodes(jump.quadrature)
    β = QuadratureRules.weights(jump.quadrature)

    μ⁻ = zeros(T, F)
    μ⁺ = zeros(T, F)
    α⁻ = zeros(T, F)
    α⁺ = zeros(T, F)

    for i in eachindex(γ)
        μ⁻[i] = evaluate_l(jump.path, γ[i])
        μ⁺[i] = evaluate_r(jump.path, γ[i])
        α⁻[i] = derivative_l(jump.path, γ[i])
        α⁺[i] = derivative_r(jump.path, γ[i])
    end

    return (
        β = SVector{F, T}(β),
        γ = SVector{F, T}(γ),
        μ⁻ = SVector{F, T}(μ⁻),
        μ⁺ = SVector{F, T}(μ⁺),
        α⁻ = SVector{F, T}(α⁻),
        α⁺ = SVector{F, T}(α⁺),
        ρ⁻ = evaluate_l(jump.path, one(T) / 2),
        ρ⁺ = evaluate_r(jump.path, one(T) / 2)
    )
end

# empty jump block for the variants that do not integrate along a path
function dgvi_jump_coefficients(::Nothing, ::Type{T}) where {T}
    (
        β = SVector{0, T}(),
        γ = SVector{0, T}(),
        μ⁻ = SVector{0, T}(),
        μ⁺ = SVector{0, T}(),
        α⁻ = SVector{0, T}(),
        α⁺ = SVector{0, T}(),
        ρ⁻ = one(T) / 2,
        ρ⁺ = one(T) / 2
    )
end

"""
Per-variant precondition on the `(basis, quadrature)` pair, checked by the generated
constructors below. The default imposes nothing; [`DGVIP0`](@ref) overrides it.
"""
check_basis_quadrature(::Type, ::Basis, ::QuadratureRule) = nothing

# The five method types share one field block, of which `DGVIPI`'s use is the widest.
# They are generated here rather than written out five times; the per-variant files add
# only the constructor, `jump!`, `residual!` and the state hooks.
for name in (:DGVI, :DGVIEXP, :DGVIP0, :DGVIP1, :DGVIPI)
    @eval begin
        struct $name{T, S, R, F, SR, BT <: Basis{T}, JT} <: DGVIMethod{T, S, R, F}
            basis::BT
            quadrature::QuadratureRule{T, R}
            jump::JT

            b::SVector{R, T}
            c::SVector{R, T}
            x::SVector{S, T}
            m::SMatrix{R, S, T, SR}
            a::SMatrix{R, S, T, SR}
            r⁻::SVector{S, T}
            r⁺::SVector{S, T}

            β::SVector{F, T}
            γ::SVector{F, T}
            μ⁻::SVector{F, T}
            μ⁺::SVector{F, T}
            α⁻::SVector{F, T}
            α⁺::SVector{F, T}
            ρ⁻::T
            ρ⁺::T

            function $name(basis::Basis{T}, quadrature::QuadratureRule{T}, jump = nothing) where {T}
                check_basis_quadrature($name, basis, quadrature)

                co = dgvi_coefficients(basis, quadrature)
                jc = jump === nothing ? dgvi_jump_coefficients(nothing, T) :
                     dgvi_jump_coefficients(jump)

                S = CompactBasisFunctions.nbasis(basis)
                R = QuadratureRules.nnodes(quadrature)
                F = length(jc.β)

                new{T, S, R, F, S * R, typeof(basis), typeof(jump)}(
                    basis, quadrature, jump,
                    co.b, co.c, co.x, co.m, co.a, co.r⁻, co.r⁺,
                    jc.β, jc.γ, jc.μ⁻, jc.μ⁺, jc.α⁻, jc.α⁺, jc.ρ⁻, jc.ρ⁺)
            end
        end
    end
end

function Base.show(io::IO, method::DGVIMethod)
    print(io, "\n", description(method), "\n")
    print(io, "  ", "="^length(description(method)), "\n\n")
    print(io, "    c  = ", method.c, "\n")
    print(io, "    b  = ", method.b, "\n")
    print(io, "    x  = ", method.x, "\n")
    print(io, "    m  = ", method.m, "\n")
    print(io, "    a  = ", method.a, "\n")
    print(io, "    r⁻ = ", method.r⁻, "\n")
    print(io, "    r⁺ = ", method.r⁺, "\n")
    if njump(method) > 0
        print(io, "    β  = ", method.β, "\n")
        print(io, "    γ  = ", method.γ, "\n")
        print(io, "    μ⁻ = ", method.μ⁻, "\n")
        print(io, "    μ⁺ = ", method.μ⁺, "\n")
        print(io, "    α⁻ = ", method.α⁻, "\n")
        print(io, "    α⁺ = ", method.α⁺, "\n")
        print(io, "    ρ⁻ = ", method.ρ⁻, "\n")
        print(io, "    ρ⁺ = ", method.ρ⁺, "\n")
    end
    print(io, "\n")
end

function Base.show(io::IO, int::GeometricIntegrator{<:DGVIMethod})
    print(io, "\n", description(method(int)), " with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, string(method(int)))
end

@doc raw"""
State that a DGVI variant carries from one step to the next.

`DGVI` is a genuine ``(q_n, p_n)`` map and uses none of it; the other four propagate
one or two of the jump values, which are *not* part of the `SolutionStep`. This lives
in the `DT` cache only — it is constant with respect to the nonlinear solver unknowns
and must never be dual-typed — and is therefore always reached through
[`dgvi_state`](@ref) rather than by naming the field.

The state is seeded once, on the first step, from the initial condition. A DGVI
`GeometricIntegrator` is therefore single-run: restarting one from different initial
conditions would carry over the jump values of the previous run. Build a fresh integrator
instead.
"""
mutable struct DGVIState{DT}
    initialised::Bool
    q⁻::Vector{DT}
    q⁺::Vector{DT}
    θ⁻::Vector{DT}

    DGVIState{DT}(D) where {DT} = new(false, zeros(DT, D), zeros(DT, D), zeros(DT, D))
end

@doc raw"""
Cache shared by all five DGVI variants.

The field set is the union of what the variants actually *read*; the legacy
per-variant caches allocated roughly twice as many vectors, most of which were
computed and never used. `Φ`, `Λ`, `Θ`, `Θ̄`, `G` and `Ḡ` hold the per-node values of
the jump path integral and are empty (`F = 0`) for every variant but
[`DGVIPI`](@ref).
"""
struct DGVICache{ST, D, S, R, F, N} <: IODEIntegratorCache{ST}
    x::Vector{ST}

    # degrees of freedom and stage values
    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    # initial-guess temporaries
    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}

    # permanent zeros, for the unused velocity slot of `g`
    z::Vector{ST}

    # boundary values
    q⁺::Vector{ST}
    q̄::Vector{ST}
    q̄⁻::Vector{ST}
    q̄⁺::Vector{ST}

    # averages and jumps
    ϕ::Vector{ST}
    ϕ̅::Vector{ST}
    λ::Vector{ST}
    λ⁺::Vector{ST}
    λ̄::Vector{ST}
    λ̄⁻::Vector{ST}

    # one-forms
    θ::Vector{ST}
    θ⁺::Vector{ST}
    Θ̅::Vector{ST}
    Θ̅⁻::Vector{ST}
    Θ̅⁺::Vector{ST}

    # projections ∇ϑᵀ·λ
    g::Vector{ST}
    g⁺::Vector{ST}
    ḡ::Vector{ST}
    ḡ⁻::Vector{ST}
    h⁺::Vector{ST}
    h̅⁻::Vector{ST}

    # momentum output
    p̄::Vector{ST}

    # jump path integral (empty unless F > 0)
    Φ::Vector{Vector{ST}}
    Φ̄::Vector{Vector{ST}}
    Λ::Vector{Vector{ST}}
    Λ̄::Vector{Vector{ST}}
    Θ::Vector{Vector{ST}}
    Θ̄::Vector{Vector{ST}}
    G::Vector{Vector{ST}}
    Ḡ::Vector{Vector{ST}}

    state::DGVIState{ST}

    # `new` takes 42 positional arguments, 27 of them interchangeable `Vector{ST}`s of
    # length D, so a misordering would not be caught by the compiler. They are therefore
    # grouped and labelled in field order, so that inserting or removing one shows up as a
    # reviewable diff rather than silently shifting everything after it.
    function DGVICache{ST, D, S, R, F, N}() where {ST, D, S, R, F, N}
        v() = zeros(ST, D)
        stage(n) = create_internal_stage_vector(ST, D, n)

        new(zeros(ST, N),                                   # x
            stage(S), stage(R), stage(R), stage(R), stage(R),  # X Q V P F
            v(), v(), v(), v(),                             # q̃ p̃ ṽ f̃
            v(),                                            # z
            v(), v(), v(), v(),                             # q⁺ q̄ q̄⁻ q̄⁺
            v(), v(), v(), v(), v(), v(),                    # ϕ ϕ̅ λ λ⁺ λ̄ λ̄⁻
            v(), v(), v(), v(), v(),                         # θ θ⁺ Θ̅ Θ̅⁻ Θ̅⁺
            v(), v(), v(), v(), v(), v(),                    # g g⁺ ḡ ḡ⁻ h⁺ h̅⁻
            v(),                                            # p̄
            stage(F), stage(F), stage(F), stage(F),           # Φ Φ̄ Λ Λ̄
            stage(F), stage(F), stage(F), stage(F),           # Θ Θ̄ G Ḡ
            DGVIState{ST}(D))
    end
end

nlsolution(cache::DGVICache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::DGVIMethod; kwargs...) where {ST}
    D = length(vec(initial_conditions(problem).q))
    DGVICache{
        ST, D, nbasis(method), nnodes(method), njump(method), solversize(method, problem)}(;
        kwargs...)
end

@inline function CacheType(ST, problem::AbstractProblemIODE, method::DGVIMethod)
    D = length(vec(initial_conditions(problem).q))
    DGVICache{
        ST, D, nbasis(method), nnodes(method), njump(method), solversize(method, problem)}
end

function internal_variables(method::DGVIMethod, problem::AbstractProblemIODE{
        DT, TT}) where {DT, TT}
    D = length(vec(initial_conditions(problem).q))
    R = nnodes(method)

    Q = create_internal_stage_vector(DT, D, R)
    V = create_internal_stage_vector(DT, D, R)
    P = create_internal_stage_vector(DT, D, R)

    (Q = Q, V = V, P = P, q⁻ = zeros(DT, D), q⁺ = zeros(DT, D))
end

# Names `cache.state` directly rather than going through `dgvi_state`, because it receives
# a cache and not an integrator. Safe: the framework always calls this with the `DT` cache
# (`copy_internal_variables!(solstep, cache(int))` in GeometricIntegratorsBase).
function copy_internal_variables!(solstep::SolutionStep, cache::DGVICache)
    haskey(internal(solstep), :Q) && copyto!(internal(solstep).Q, cache.Q)
    haskey(internal(solstep), :V) && copyto!(internal(solstep).V, cache.V)
    haskey(internal(solstep), :P) && copyto!(internal(solstep).P, cache.P)
    haskey(internal(solstep), :q⁻) && copyto!(internal(solstep).q⁻, cache.state.q⁻)
    haskey(internal(solstep), :q⁺) && copyto!(internal(solstep).q⁺, cache.state.q⁺)
end

@doc raw"""
The jump values a DGVI variant carries from one step to the next, as a [`DGVIState`](@ref).

Always reach the state through this function. `DGVIState` is a field of the `ST`-parameterised
[`DGVICache`](@ref), so `cache(int, ST).state` also exists — but as a *separate*,
dual-typed object, so writing to it would be silently lost. Going through `dgvi_state`
keeps that impossible to get wrong.
"""
@inline dgvi_state(int::GeometricIntegrator{<:DGVIMethod}) = cache(int).state

"""
Seed the carried-over jump values on the first step. The generic version starts from a
continuous trajectory, `qₙ⁻ = qₙ⁺ = qₙ`, and `ϑ(qₙ⁻) = pₙ`; variants override it.
"""
function initialise_state!(sol, params, int::GeometricIntegrator{<:DGVIMethod})
    local state = dgvi_state(int)
    state.initialised && return
    state.q⁻ .= sol.q
    state.q⁺ .= sol.q
    state.θ⁻ .= sol.p
    state.initialised = true
    return
end

"""
Carry the jump values of this step over to the next one. Overridden per variant.
"""
update_state!(int::GeometricIntegrator{<:DGVIMethod}, DT) = nothing

function initial_guess!(sol, history, params, int::GeometricIntegrator{<:DGVIMethod})
    local x = nlsolution(int)
    local D = length(cache(int).q̃)
    local S = nbasis(method(int))

    # TODO: this initialises the degrees of freedom with the solution itself; an L2
    # projection of q onto the basis would be more accurate.
    for i in eachindex(basis(method(int)))
        soltmp = (
            t = sol.t + timestep(int) * (method(int).x[i] - 1),
            q = cache(int).q̃,
            p = cache(int).p̃,
            q̇ = cache(int).ṽ,
            ṗ = cache(int).f̃
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:D
            x[D * (i - 1) + k] = cache(int).q̃[k]
        end
    end

    # the trailing block(s) are seeded with the solution at the end of the step
    soltmp = (
        t = sol.t,
        q = cache(int).q̃,
        p = cache(int).p̃,
        q̇ = cache(int).ṽ,
        ṗ = cache(int).f̃
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    for i in (S + 1):(length(x) ÷ D)
        for k in 1:D
            x[D * (i - 1) + k] = cache(int).q̃[k]
        end
    end

    initialise_state!(sol, params, int)
end

"""
Compute the solution at the quadrature nodes and the two boundary reconstructions
`qₙ⁺ = r⁺·X` and `qₙ₊₁⁻ = r⁻·X`.
"""
function components_q!(sol, params, int::GeometricIntegrator{<:DGVIMethod}, ST)
    local C = cache(int, ST)
    local M = method(int)

    for i in eachindex(C.Q)
        for k in eachindex(C.Q[i])
            y = zero(ST)
            for j in eachindex(C.X)
                y += M.m[i, j] * C.X[j][k]
            end
            C.Q[i][k] = y
        end
    end

    for k in eachindex(C.q⁺, C.q̄⁻)
        y⁺ = y⁻ = zero(ST)
        for i in eachindex(C.X)
            y⁺ += M.r⁺[i] * C.X[i][k]
            y⁻ += M.r⁻[i] * C.X[i][k]
        end
        C.q⁺[k] = y⁺
        C.q̄⁻[k] = y⁻
    end
end

"Compute the velocities at the quadrature nodes, `V = a X / h`."
function components_v!(sol, params, int::GeometricIntegrator{<:DGVIMethod}, ST)
    local C = cache(int, ST)
    local M = method(int)

    for i in eachindex(C.V)
        for k in eachindex(C.V[i])
            y = zero(ST)
            for j in eachindex(C.X)
                y += M.a[i, j] * C.X[j][k]
            end
            C.V[i][k] = y / timestep(int)
        end
    end
end

"Compute the one-form and the forces at the quadrature nodes."
function components_p!(sol, params, int::GeometricIntegrator{<:DGVIMethod}, ST)
    local C = cache(int, ST)

    for i in eachindex(C.P, C.F)
        tᵢ = sol.t + timestep(int) * (method(int).c[i] - 1)
        equations(int).ϑ(C.P[i], tᵢ, C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], tᵢ, C.Q[i], C.V[i], params)
    end
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:DGVIMethod}) where {ST}
    local C = cache(int, ST)
    local D = length(C.q̃)
    local S = nbasis(method(int))

    # copy x to the degrees of freedom
    for i in eachindex(C.X)
        for k in eachindex(C.X[i])
            C.X[i][k] = x[D * (i - 1) + k]
        end
    end

    # the first trailing block is always qₙ₊₁
    for k in eachindex(C.q̄)
        C.q̄[k] = x[D * S + k]
    end

    # the second trailing block, where present, is qₙ₊₁⁺
    if nclosure(method(int)) > 1
        for k in eachindex(C.q̄⁺)
            C.q̄⁺[k] = x[D * (S + 1) + k]
        end
    end

    components_q!(sol, params, int, ST)
    components_v!(sol, params, int, ST)
    components_p!(sol, params, int, ST)

    jump!(sol, params, int, ST)
end

"""
Interior contribution to the discrete Euler-Lagrange equations, shared by all variants:
`b[i] = Σⱼ bⱼ (h mⱼᵢ Fⱼ + aⱼᵢ Pⱼ)`. The variant's `residual!` adds its flux terms on
top and writes the trailing block(s).
"""
function residual_interior!(b::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:DGVIMethod}) where {ST}
    local C = cache(int, ST)
    local M = method(int)
    local D = length(C.q̃)

    for i in eachindex(M.r⁻, M.r⁺)
        for k in 1:D
            z = zero(ST)
            for j in eachindex(C.P, C.F)
                z += M.b[j] * M.m[j, i] * C.F[j][k] * timestep(int)
                z += M.b[j] * M.a[j, i] * C.P[j][k]
            end
            b[D * (i - 1) + k] = z
        end
    end
end

function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:DGVIMethod}) where {ST}
    @assert axes(x) == axes(b)
    components!(x, sol, params, int)
    residual!(b, sol, params, int)
end

function update!(sol, params, int::GeometricIntegrator{<:DGVIMethod}, DT)
    sol.q .= cache(int, DT).q̄
    sol.p .= cache(int, DT).p̄
    update_state!(int, DT)
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:DGVIMethod}) where {DT}
    components!(x, sol, params, int)
    update!(sol, params, int, DT)
end

function integrate_step!(sol, history, params, int::GeometricIntegrator{
        <:DGVIMethod, <:AbstractProblemIODE})
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (
        sol, params, int))
    check_solver_status(solverstatus, int)
    update!(sol, params, nlsolution(int), int)
end

@doc raw"""
Abstract supertype of the Continuous Galerkin Variational Integrators.

A CGVI discretises the action of a Lagrangian by a *continuous* piecewise-polynomial
trajectory. Within each interval,
```math
q_h(t) \vert_{[t_{n}, t_{n+1}]} = \sum \limits_{j=1}^{S} X_{n,j} \, \varphi_{j} \bigg( \frac{t - t_{n}}{t_{n+1} - t_{n}} \bigg) ,
```
and with a quadrature rule of ``R`` nodes ``c_i`` and weights ``b_i`` the mass and
derivative matrices and the boundary reconstruction vectors are
```math
m_{ij} = \varphi_j (c_i) ,
\qquad
a_{ij} = \varphi_j' (c_i) ,
\qquad
r^{0}_{j} = \varphi_j (0) ,
\qquad
r^{1}_{j} = \varphi_j (1) ,
```
so that the stage values and velocities are
```math
Q_{n,i} = m_{ij} X_{n,j} ,
\qquad
V_{n,i} = \frac{a_{ij}}{h} X_{n,j} ,
```
with ``P_{n,i} = \vartheta (Q_{n,i}, V_{n,i})`` and ``F_{n,i} = f (Q_{n,i}, V_{n,i})``.
All variants share this construction and the interior contribution to the discrete
Euler-Lagrange equations,
```math
\sum \limits_{i=1}^{R} b_i \big[ h \, m_{ij} \, F_{n,i} + a_{ij} \, P_{n,i} \big] ,
```
and differ only in how continuity across the interval boundaries is imposed, and hence
in which quantities are unknowns of the nonlinear system.

The two variants are [`CGVI`](@ref), which adds the continuity constraint to the action
through Lagrange multipliers, and [`CGVINodal`](@ref), which builds it into the basis.
See `docs/src/integrators/cgvi.md`.
"""
abstract type CGVIMethod{T, S, R} <: LODEMethod end

isexplicit(::Union{CGVIMethod, Type{<:CGVIMethod}}) = false
isimplicit(::Union{CGVIMethod, Type{<:CGVIMethod}}) = true
issymmetric(::Union{CGVIMethod, Type{<:CGVIMethod}}) = missing
issymplectic(::Union{CGVIMethod, Type{<:CGVIMethod}}) = true

isiodemethod(::Union{CGVIMethod, Type{<:CGVIMethod}}) = true

default_solver(::CGVIMethod) = Newton()
default_iguess(::CGVIMethod) = HermiteExtrapolation()

basis(method::CGVIMethod) = method.basis
quadrature(method::CGVIMethod) = method.quadrature

nbasis(::CGVIMethod{T, S, R}) where {T, S, R} = S
nnodes(::CGVIMethod{T, S, R}) where {T, S, R} = R

"""
Compute the shared CGVI coefficient block `(b, c, x, m, a, r₀, r₁)` from a basis and a
quadrature rule.
"""
function cgvi_coefficients(basis::Basis{T}, quadrature::QuadratureRule{T}) where {T}
    R = QuadratureRules.nnodes(quadrature)
    S = CompactBasisFunctions.nbasis(basis)

    quad_nodes = QuadratureRules.nodes(quadrature)
    quad_weights = QuadratureRules.weights(quadrature)

    m = zeros(T, R, S)
    a = zeros(T, R, S)
    r₀ = zeros(T, S)
    r₁ = zeros(T, S)

    for i in eachindex(basis)
        r₀[i] = basis[zero(T), i]
        r₁[i] = basis[one(T), i]
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
        r₀ = SVector{S, T}(r₀),
        r₁ = SVector{S, T}(r₁)
    )
end

function Base.show(io::IO, method::CGVIMethod)
    local title = description(method)
    print(io, "\n")
    print(io, "  ", title, "\n")
    print(io, "  ", "="^length(title), "\n")
    print(io, "\n")
    print(io, "    c  = ", method.c, "\n")
    print(io, "    b  = ", method.b, "\n")
    print(io, "    x  = ", method.x, "\n")
    print(io, "    m  = ", method.m, "\n")
    print(io, "    a  = ", method.a, "\n")
    print(io, "    r₀ = ", method.r₀, "\n")
    print(io, "    r₁ = ", method.r₁, "\n")
    print(io, "\n")
end

@doc raw"""
Cache shared by both CGVI variants.

`X` holds the `S` basis coefficients and `Q`, `V`, `P`, `F` the stage values at the `R`
quadrature nodes; the variants agree on all of these and differ only in the length of the
nonlinear solution vector `x`, which comes from `solversize`. `q̃`, `p̃`, `ṽ`, `f̃` are the
initial-guess temporaries; [`CGVI`](@ref) additionally uses `q̃` and `p̃` to carry the
reconstructed endpoint values into `update!`.

The number of degrees of freedom `D` and the solver size are constructor *arguments*, not
type parameters, so that `CacheType` below stays a function of the method alone. It is
evaluated on every `cache(int, ST)` — `CacheDict` type-asserts its `getindex` on it — and
a `CacheType` that reads a value off the problem does not constant-fold, which leaves the
cache inferred abstractly and makes the Newton hot path box on every stage access.
"""
struct CGVICache{ST, S, R} <: IODEIntegratorCache{ST}
    x::Vector{ST}

    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    F::Vector{Vector{ST}}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}

    function CGVICache{ST, S, R}(D::Int, N::Int) where {ST, S, R}
        v() = zeros(ST, D)
        stage(n) = create_internal_stage_vector(ST, D, n)

        new(zeros(ST, N),
            stage(S), stage(R), stage(R), stage(R), stage(R),
            v(), v(), v(), v())
    end
end

nlsolution(cache::CGVICache) = cache.x

"Number of degrees of freedom of the problem the cache was built for."
ndofs(cache::CGVICache) = length(cache.q̃)

function Cache{ST}(problem::AbstractProblemIODE, method::CGVIMethod; kwargs...) where {ST}
    D = length(vec(initial_conditions(problem).q))
    CGVICache{ST, nbasis(method), nnodes(method)}(D, solversize(method, problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblemIODE, method::CGVIMethod) = CGVICache{
    ST, nbasis(method), nnodes(method)}

"""
Evaluate the initial guess at the relative node `τ ∈ [0,1]` of the current step, into the
cache temporaries `q̃`, `p̃`, `ṽ`, `f̃`.
"""
function initial_guess_at!(sol, history, int::GeometricIntegrator{<:CGVIMethod}, τ)
    soltmp = (
        t = sol.t + timestep(int) * (τ - 1),
        q = cache(int).q̃,
        p = cache(int).p̃,
        q̇ = cache(int).ṽ,
        ṗ = cache(int).f̃
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))
end

"""
Seed consecutive position blocks of the nonlinear solution vector `x` from the initial
guess evaluated at `nodes`, writing block `i` at `x[D*(i-1)+k]`.

The variants pass different node lists — [`CGVI`](@ref) all `S` basis nodes,
[`CGVINodal`](@ref) the `S-1` nodes left free once the initial condition has pinned the
first — but the block layout is the same, and it is the layout `components!` and
`residual!` assume.
"""
function initial_guess_positions!(
        x, sol, history, int::GeometricIntegrator{<:CGVIMethod}, nodes)
    local D = ndofs(cache(int))

    # TODO: here we should not initialise with the solution q but with the degree of freedom x,
    # obtained e.g. from an L2 projection of q onto the basis

    for (i, τ) in pairs(nodes)
        initial_guess_at!(sol, history, int, τ)

        for k in 1:D
            x[D * (i - 1) + k] = cache(int).q̃[k]
        end
    end
end

"Compute the solution at the quadrature nodes, `Q = m X`."
function components_q!(int::GeometricIntegrator{<:CGVIMethod}, ::Type{ST}) where {ST}
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
end

"Compute the velocities at the quadrature nodes, `V = a X / h`."
function components_v!(int::GeometricIntegrator{<:CGVIMethod}, ::Type{ST}) where {ST}
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

"Compute the one-form `P = ϑ(Q,V)` and the forces `F = f(Q,V)` at the quadrature nodes."
function components_p!(sol, params, int::GeometricIntegrator{<:CGVIMethod}, ::Type{ST}) where {ST}
    local C = cache(int, ST)

    for i in eachindex(C.Q, C.V, C.P, C.F)
        tᵢ = sol.t + timestep(int) * (method(int).c[i] - 1)
        equations(int).ϑ(C.P[i], tᵢ, C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], tᵢ, C.Q[i], C.V[i], params)
    end
end

function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params,
        int::GeometricIntegrator{<:CGVIMethod}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:CGVIMethod}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end

function integrate_step!(sol, history, params, int::GeometricIntegrator{
        <:CGVIMethod, <:AbstractProblemIODE})
    # call nonlinear solver and act on the outcome it reports
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (
        sol, params, int))
    check_solver_status(solverstatus, int)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end

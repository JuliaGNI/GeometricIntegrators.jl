@doc raw"""
    CGVI_II(basis, quadrature)

Continuous Galerkin Variational Integrator based on a **Type II generating function**.

## Mathematical background

On each time step $[t_n, t_{n+1}]$ of length $h$, the configuration trajectory is
approximated by a polynomial expressed in a basis $\{\varphi_j\}_{j=1}^{S}$:

```math
q_h(\tau) = \sum_{j=1}^{S} X_j \, \varphi_j(\tau), \qquad \tau \in [0,1].
```

The basis must have nodes at both endpoints of $[0,1]$ (Lobatto-type), so that
$\varphi_1(0) = 1$ and $\varphi_S(1) = 1$. This ensures that:

* $X_1 = q_n$ can be pinned to the current position (left boundary condition).
* $X_S$ equals $q_{n+1}$ exactly (no endpoint reconstruction needed).

The **free degrees of freedom** are $X_2, \ldots, X_S$ (size $D \times (S-1)$).

## Discrete equations

Internal stages at quadrature nodes $\{c_i, b_i\}_{i=1}^{R}$:

```math
Q_i = \sum_j m_{ij} X_j, \qquad
V_i = \frac{1}{h}\sum_j a_{ij} X_j,
```

where $m_{ij} = \varphi_j(c_i)$ and $a_{ij} = \varphi_j'(c_i)$.

The nonlinear system has $D \times (S-1)$ equations:

* **Left boundary** ($D$ equations, determines consistency with $p_n$):
```math
p_n + \sum_j b_j \bigl[ m_{j1}\, F_j \, h + a_{j1}\, P_j \bigr] = 0
```

* **Interior Euler–Lagrange** ($(S-2) \times D$ equations, for $i = 2, \ldots, S-1$):
```math
\sum_j b_j \bigl[ m_{j,i}\, F_j \, h + a_{j,i}\, P_j \bigr] = 0
```

where $P_i = \vartheta(Q_i, V_i)$ and $F_i = f(Q_i, V_i)$.

## Update step

After solving the nonlinear system:

```math
q_{n+1} = X_S, \qquad
p_{n+1} = \sum_j b_j \bigl[ m_{jS}\, F_j \, h + a_{jS}\, P_j \bigr].
```

The momentum $p_{n+1}$ is obtained as the discrete right-endpoint Legendre transform.

## Fields

* `b`: quadrature weights $b_i$
* `c`: quadrature nodes $c_i$
* `x`: basis nodes (must include 0 and 1, e.g. Lobatto points)
* `m`: mass matrix, $m_{ij} = \varphi_j(c_i)$
* `a`: derivative matrix, $a_{ij} = \varphi_j'(c_i)$
* `r₀`: reconstruction coefficients at $\tau=0$, $r_{0j} = \varphi_j(0)$
* `r₁`: reconstruction coefficients at $\tau=1$, $r_{1j} = \varphi_j(1)$

## Requirements

The basis nodes must include the interval endpoints (e.g. Lagrange basis on
Lobatto–Legendre nodes). Using interior-only nodes (e.g. Gauss–Legendre) will
produce incorrect results because $X_1 \neq q(0)$ and $X_S \neq q(1)$.

## Example

```julia
Q = LobattoLegendreQuadrature(4)
B = Lagrange(QuadratureRules.nodes(Q))
method = CGVI_II(B, Q)
sol = integrate(problem, method)
```

"""
struct CGVI_II{T,NBASIS,NNODES,NDOF,basisType<:Basis{T}} <: LODEMethod
    basis::basisType
    quadrature::QuadratureRule{T,NNODES}

    b::SVector{NNODES,T}
    c::SVector{NNODES,T}

    x::SVector{NBASIS,T}

    m::SMatrix{NNODES,NBASIS,T,NDOF}
    a::SMatrix{NNODES,NBASIS,T,NDOF}

    r₀::SVector{NBASIS,T}
    r₁::SVector{NBASIS,T}

    function CGVI_II(basis::Basis{T}, quadrature::QuadratureRule{T}) where {T}
        # get number of quadrature nodes and number of basis functions
        NNODES = QuadratureRules.nnodes(quadrature)
        NBASIS = CompactBasisFunctions.nbasis(basis)

        # get quadrature nodes and weights
        quad_weights = QuadratureRules.weights(quadrature)
        quad_nodes = QuadratureRules.nodes(quadrature)

        # compute coefficients
        r₀ = zeros(T, NBASIS)
        r₁ = zeros(T, NBASIS)
        m = zeros(T, NNODES, NBASIS)
        a = zeros(T, NNODES, NBASIS)

        for i in eachindex(basis)
            r₀[i] = basis[zero(T), i]
            r₁[i] = basis[one(T), i]
            for j in eachindex(quad_nodes)
                m[j, i] = basis[quad_nodes[j], i]
                a[j, i] = basis'[quad_nodes[j], i]
            end
        end

        new{T,NBASIS,NNODES,NBASIS * NNODES,typeof(basis)}(basis, quadrature, quad_weights, quad_nodes, CompactBasisFunctions.grid(basis), m, a, r₀, r₁)
    end
end

basis(method::CGVI_II) = method.basis
quadrature(method::CGVI_II) = method.quadrature

nbasis(::CGVI_II{T,NB,NN}) where {T,NB,NN} = NB
nnodes(::CGVI_II{T,NB,NN}) where {T,NB,NN} = NN

isexplicit(::Union{CGVI_II,Type{<:CGVI_II}}) = false
isimplicit(::Union{CGVI_II,Type{<:CGVI_II}}) = true
issymmetric(::Union{CGVI_II,Type{<:CGVI_II}}) = missing
issymplectic(::Union{CGVI_II,Type{<:CGVI_II}}) = true

isiodemethod(::Union{CGVI_II,Type{<:CGVI_II}}) = true

default_solver(::CGVI_II) = Newton()
default_iguess(::CGVI_II) = HermiteExtrapolation()

function Base.show(io::IO, method::CGVI_II)
    print(io, "\n")
    print(io, "  Continuous Galerkin Variational Integrator", "\n")
    print(io, "  ==========================================", "\n")
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


struct CGVI_IICache{ST,S,R} <: IODEIntegratorCache{ST}
    x::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    s̃::Vector{ST}

    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}


    function CGVI_IICache{ST,S,R}(ics) where {ST,S,R}
        D = length(vec(ics.q))
        x = zeros(ST, D * (S-1))

        # create temporary vectors
        q̃ = zeros(ST, D)
        p̃ = zeros(ST, D)
        ṽ = zeros(ST, D)
        f̃ = zeros(ST, D)
        s̃ = zeros(ST, D)

        # create internal stage vectors
        X = create_internal_stage_vector(ST, D, S)
        Q = create_internal_stage_vector(ST, D, R)
        P = create_internal_stage_vector(ST, D, R)
        V = create_internal_stage_vector(ST, D, R)
        F = create_internal_stage_vector(ST, D, R)

        new(x, q̃, p̃, ṽ, f̃, s̃, X, Q, P, V, F)
    end
end

nlsolution(cache::CGVI_IICache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::CGVI_II; kwargs...) where {ST}
    CGVI_IICache{ST,nbasis(method),nnodes(method)}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblemIODE, method::CGVI_II) = CGVI_IICache{ST,nbasis(method),nnodes(method)}


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CGVI_II})
    # set some local variables for convenience
    local D = length(cache(int).q̃)
    local S = nbasis(method(int))
    local x = nlsolution(int)

    # TODO: here we should not initialise with the solution q but with the degree of freedom x,
    # obtained e.g. from an L2 projection of q onto the basis

    for i in 1:length(method(int).x)-1
        soltmp = (
            t=sol.t + timestep(int) * (method(int).x[i+1] - 1),
            q=cache(int).q̃,
            p=cache(int).p̃,
            q̇=cache(int).ṽ,
            ṗ=cache(int).f̃,
        )
        solutionstep!(soltmp, history, problem(int), iguess(int))

        for k in 1:D
            x[D*(i-1)+k] = cache(int).q̃[k]
        end
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CGVI_II}) where {ST}
    # set some local variables for convenience and clarity
    local C = cache(int, ST)
    local D = length(C.q̃)
    local S = nbasis(method(int))
    local q̄ = sol.q

    for d in 1:D
        C.X[1][d] = q̄[d]
    end

    # copy x to X
    for s in 1:S-1
        for d in 1:D
            C.X[s+1][d] = x[D*(s-1)+d]
        end
    end

    # compute Q
    for i in eachindex(C.Q)
        for k in eachindex(C.Q[i])
            y = zero(ST)
            for j in eachindex(C.X)
                y += method(int).m[i, j] * C.X[j][k]
            end
            C.Q[i][k] = y
        end
    end

    # compute V
    for i in eachindex(C.V)
        for k in eachindex(C.V[i])
            y = zero(ST)
            for j in eachindex(C.X)
                y += method(int).a[i, j] * C.X[j][k]
            end
            C.V[i][k] = y / timestep(int)
        end
    end

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in eachindex(C.Q, C.V, C.P, C.F)
        tᵢ = sol.t + timestep(int) * (method(int).c[i] - 1)
        equations(int).ϑ(C.P[i], tᵢ, C.Q[i], C.V[i], params)
        equations(int).f(C.F[i], tᵢ, C.Q[i], C.V[i], params)
    end
end


function residual!(b::Vector{ST}, sol, params, int::GeometricIntegrator{<:CGVI_II}) where {ST}
    # set some local variables for convenience and clarity
    local C = cache(int, ST)
    local D = length(C.q̃)
    local S = nbasis(method(int))
    local p̄ = sol.p

    for k in eachindex(p̄)
        z = zero(ST)
        for j in eachindex(C.P, C.F)
            z += method(int).b[j] * C.F[j][k] * method(int).m[j, 1] * timestep(int)
            z += method(int).b[j] * C.P[j][k] * method(int).a[j, 1]
        end
        b[k] = p̄[k] + z
    end

    # compute b = - [(P-AF)]
    for i in 1:S-2  
        for k in 1:D 
            z = zero(ST)
            for j in eachindex(C.P, C.F) # quad_nodes index 
                z += method(int).b[j] * method(int).m[j, i+1] * C.F[j][k] * timestep(int)
                z += method(int).b[j] * method(int).a[j, i+1] * C.P[j][k]
            end
            b[D + D*(i-1)+k] = z
        end
    end
end


# Compute stages of Variational Partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CGVI_II}) where {ST}
    # check that x and b are compatible
    @assert axes(x) == axes(b)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, int::GeometricIntegrator{<:CGVI_II}, DT)
    local C = cache(int, DT)
    local D = length(C.q̃)
    local S = nbasis(method(int))
    local h = timestep(int)

    sol.q .= C.X[end]

    for k in 1:D
        z = zero(DT)
        for j in 1:nnodes(method(int))
            z += method(int).b[j] * C.F[j][k] * method(int).m[j, S] * h
            z += method(int).b[j] * C.P[j][k] * method(int).a[j, S]
        end
        sol.p[k] = z
    end
end

function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:CGVI_II}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    update!(sol, params, int, DT)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:CGVI_II,<:AbstractProblemIODE})
    # call nonlinear solver
    # solve!(nlsolution(int), (b, x) -> residual!(b, x, sol, params, int), solver(int))
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
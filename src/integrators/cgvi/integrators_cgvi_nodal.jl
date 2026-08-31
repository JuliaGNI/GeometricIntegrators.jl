@doc raw"""
Continuous Galerkin Variational Integrator on a *nodal* basis.

In the Galerkin framework of [OberBloebaum:2015](@cite) the continuity of the trajectory at
the interval boundaries is enforced by the approximation polynomial itself rather than by a
constraint added to the action: the basis is interpolatory with nodes
``0 = x_1 < x_2 < \dots < x_S = 1``, so that ``\varphi_1(0) = 1`` and ``\varphi_S(1) = 1``
and the coefficients *are* nodal values. The first is therefore pinned to the known ``q_n``,
the last *is* ``q_{n+1}``, and the momentum is computed explicitly rather than solved for.

That leaves `D*(S-1)` unknowns in the nonlinear system, against the `D*(S+1)` of
[`CGVI`](@ref), which solves for all `S` coefficients plus the endpoint momentum. No
Lagrange multipliers are needed here, since the constraint they enforce in `CGVI` holds
by construction.

## Requirements

The basis must be interpolatory with nodes at both ends of the interval — a Lagrange basis
on Lobatto-Legendre nodes satisfies this. Interior-only nodes (e.g. Gauss-Legendre) would
give ``X_1 \neq q(0)`` and ``X_S \neq q(1)`` and hence silently wrong results, so the
constructor rejects them.

## Example

```julia
Q = LobattoLegendreQuadrature(4)
B = Lagrange(QuadratureRules.nodes(Q))
sol = integrate(problem, CGVINodal(B, Q))
```

## Fields

* `b`: weights of the quadrature rule
* `c`: nodes of the quadrature rule
* `x`: nodes of the basis
* `m`: mass matrix
* `a`: derivative matrix
* `r₀`: reconstruction coefficients at the beginning of the interval
* `r₁`: reconstruction coefficients at the end of the interval

"""
struct CGVINodal{T, NBASIS, NNODES, NMA, basisType <: Basis{T}} <:
       CGVIMethod{T, NBASIS, NNODES}
    basis::basisType
    quadrature::QuadratureRule{T, NNODES}

    b::SVector{NNODES, T}
    c::SVector{NNODES, T}

    x::SVector{NBASIS, T}

    m::SMatrix{NNODES, NBASIS, T, NMA}
    a::SMatrix{NNODES, NBASIS, T, NMA}

    r₀::SVector{NBASIS, T}
    r₁::SVector{NBASIS, T}

    function CGVINodal(basis::Basis{T}, quadrature::QuadratureRule{T}) where {T}
        NNODES = QuadratureRules.nnodes(quadrature)
        NBASIS = CompactBasisFunctions.nbasis(basis)

        coeffs = cgvi_coefficients(basis, quadrature)

        # The formulation reads q(0) off X[1] and q(1) off X[S] instead of reconstructing
        # them, which is only valid when φ₁(0) = φ_S(1) = 1 and every other basis function
        # vanishes there. A basis on interior nodes passes every other check and produces
        # a wrong trajectory, so it is rejected here rather than only documented.
        e₀ = SVector{NBASIS, T}(ntuple(i -> T(i == 1), NBASIS))
        e₁ = SVector{NBASIS, T}(ntuple(i -> T(i == NBASIS), NBASIS))
        @assert isapprox(coeffs.r₀, e₀; atol = sqrt(eps(T))) &&
                isapprox(coeffs.r₁, e₁; atol = sqrt(eps(T))) (
            "CGVINodal requires an interpolatory basis with nodes at both ends of the interval, " *
            "e.g. Lagrange(QuadratureRules.nodes(LobattoLegendreQuadrature(s))). " *
            "Got φ(0) = $(coeffs.r₀) and φ(1) = $(coeffs.r₁), " *
            "expected $e₀ and $e₁.")

        new{T, NBASIS, NNODES, NBASIS * NNODES, typeof(basis)}(
            basis, quadrature, coeffs.b, coeffs.c, coeffs.x,
            coeffs.m, coeffs.a, coeffs.r₀, coeffs.r₁)
    end
end

function GeometricBase.description(::CGVINodal)
    "Continuous Galerkin Variational Integrator (nodal basis)"
end

# all `S` basis coefficients but the first, which the initial condition pins to qₙ
function solversize(method::CGVINodal, problem::AbstractProblemIODE)
    length(vec(initial_conditions(problem).q)) * (nbasis(method) - 1)
end

function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CGVINodal})
    # the first coefficient is pinned to q̄ and is not an unknown, so only the remaining
    # basis nodes get a block
    initial_guess_positions!(
        nlsolution(int), sol, history, int, Base.tail(Tuple(method(int).x)))
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CGVINodal}) where {ST}
    local C = cache(int, ST)
    local D = ndofs(C)
    local S = nbasis(method(int))
    local q̄ = sol.q

    for d in 1:D
        C.X[1][d] = q̄[d]
    end

    # Copy x to X. The nonlinear solution vector holds the `S-1` free basis coefficients
    # `X[2] … X[S]` with the degree of freedom running fastest: coefficient `s+1` of degree of
    # freedom `d` sits at `x[D*(s-1)+d]`. That is the layout `initial_guess!` writes and the one
    # `residual!` assumes for `b`; all three have to agree or the Jacobian picks up a zero column.
    for s in 1:(S - 1)
        for d in 1:D
            C.X[s + 1][d] = x[D * (s - 1) + d]
        end
    end

    components_q!(int, ST)
    components_v!(int, ST)
    components_p!(sol, params, int, ST)
end

function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CGVINodal}) where {ST}
    local C = cache(int, ST)
    local D = ndofs(C)
    local M = method(int)
    local S = nbasis(M)
    local p̄ = sol.p

    # the first basis function carries the momentum matching condition at the beginning of
    # the interval, since X[1] is pinned rather than free
    for k in 1:D
        z = zero(ST)
        for j in eachindex(C.P, C.F)
            z += M.b[j] * C.F[j][k] * M.m[j, 1] * timestep(int)
            z += M.b[j] * C.P[j][k] * M.a[j, 1]
        end
        b[k] = p̄[k] + z
    end

    # the interior basis functions carry the discrete Euler-Lagrange equations
    for i in 1:(S - 2)
        for k in 1:D
            z = zero(ST)
            for j in eachindex(C.P, C.F)
                z += M.b[j] * M.m[j, i + 1] * C.F[j][k] * timestep(int)
                z += M.b[j] * M.a[j, i + 1] * C.P[j][k]
            end
            b[D + D * (i - 1) + k] = z
        end
    end
end

function update!(sol, params, int::GeometricIntegrator{<:CGVINodal}, DT)
    local C = cache(int, DT)
    local D = ndofs(C)
    local M = method(int)
    local S = nbasis(M)
    local h = timestep(int)

    # The basis is interpolatory with a node at the end of the interval, so the last coefficient
    # *is* the new position — one value per degree of freedom. Read it off `C.X`, which
    # `components!` has just filled, rather than indexing into the flat solution vector.
    for k in 1:D
        sol.q[k] = C.X[S][k]
    end

    for k in 1:D
        z = zero(DT)
        for j in eachindex(C.P, C.F)
            z += M.b[j] * C.F[j][k] * M.m[j, S] * h
            z += M.b[j] * C.P[j][k] * M.a[j, S]
        end
        sol.p[k] = z
    end
end

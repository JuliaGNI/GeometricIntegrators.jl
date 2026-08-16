@doc raw"""
Continuous Galerkin Variational Integrator.

The continuity of the trajectory across the interval boundaries is imposed *weakly*, by
adding the constraints
```math
q_{h} \vert_{[t_{n}, t_{n+1}]} (t_{n}) = q_{n}
\qquad \text{and} \qquad
q_{h} \vert_{[t_{n}, t_{n+1}]} (t_{n+1}) = q_{n+1}
```
to the discrete action through Lagrange multipliers. The unknowns of the nonlinear system
are therefore all `S` basis coefficients together with the momentum at the end of the
interval, `D*(S+1)` in total, and the endpoint position is reconstructed as `q = r₁·X`.
Any basis may be used; the nodes need not include the interval endpoints.

See [`CGVINodal`](@ref) for the variant that builds continuity into the basis instead, and
[`CGVIMethod`](@ref) for the construction both share.

* `b`: weights of the quadrature rule
* `c`: nodes of the quadrature rule
* `x`: nodes of the basis
* `m`: mass matrix
* `a`: derivative matrix
* `r₀`: reconstruction coefficients at the beginning of the interval
* `r₁`: reconstruction coefficients at the end of the interval

"""
struct CGVI{T,NBASIS,NNODES,NMA,basisType<:Basis{T}} <: CGVIMethod{T,NBASIS,NNODES}
    basis::basisType
    quadrature::QuadratureRule{T,NNODES}

    b::SVector{NNODES,T}
    c::SVector{NNODES,T}

    x::SVector{NBASIS,T}

    m::SMatrix{NNODES,NBASIS,T,NMA}
    a::SMatrix{NNODES,NBASIS,T,NMA}

    r₀::SVector{NBASIS,T}
    r₁::SVector{NBASIS,T}

    function CGVI(basis::Basis{T}, quadrature::QuadratureRule{T}) where {T}
        NNODES = QuadratureRules.nnodes(quadrature)
        NBASIS = CompactBasisFunctions.nbasis(basis)

        coeffs = cgvi_coefficients(basis, quadrature)

        new{T,NBASIS,NNODES,NBASIS * NNODES,typeof(basis)}(
            basis, quadrature, coeffs.b, coeffs.c, coeffs.x, coeffs.m, coeffs.a, coeffs.r₀, coeffs.r₁)
    end
end

GeometricBase.description(::CGVI) = "Continuous Galerkin Variational Integrator"

# all `S` basis coefficients, plus the momentum at the end of the interval
solversize(method::CGVI, problem::AbstractProblemIODE) =
    length(vec(initial_conditions(problem).q)) * (nbasis(method) + 1)


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CGVI})
    local D = ndofs(cache(int))
    local S = nbasis(method(int))
    local x = nlsolution(int)

    # one position block per basis node
    initial_guess_positions!(x, sol, history, int, method(int).x)

    # the trailing block is the momentum at the end of the step
    initial_guess_at!(sol, history, int, one(eltype(method(int).x)))

    for k in 1:D
        x[D*S+k] = cache(int).p̃[k]
    end
end


function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CGVI}) where {ST}
    local C = cache(int, ST)
    local D = ndofs(C)
    local S = nbasis(method(int))

    # copy x to X
    for i in eachindex(C.X)
        for k in eachindex(C.X[i])
            C.X[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to p
    for k in eachindex(C.p̃)
        C.p̃[k] = x[D*S+k]
    end

    components_q!(int, ST)

    # reconstruct the position at the end of the interval, q = r₁·X
    for k in eachindex(C.q̃)
        y = zero(ST)
        for i in eachindex(C.X)
            y += method(int).r₁[i] * C.X[i][k]
        end
        C.q̃[k] = y
    end

    components_v!(int, ST)
    components_p!(sol, params, int, ST)
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CGVI}) where {ST}
    local C = cache(int, ST)
    local D = ndofs(C)
    local M = method(int)
    local S = nbasis(M)

    # compute b = - [(P-AF)]
    for i in eachindex(M.r₀, M.r₁)
        for k in eachindex(C.p̃)#, sol.p # TODO
            z = zero(ST)
            for j in eachindex(C.P, C.F)
                z += M.b[j] * M.m[j, i] * C.F[j][k] * timestep(int)
                z += M.b[j] * M.a[j, i] * C.P[j][k]
            end
            b[D*(i-1)+k] = (M.r₁[i] * C.p̃[k] - M.r₀[i] * sol.p[k]) - z
        end
    end

    # compute b = - [(q-r₀Q)]
    for k in eachindex(sol.q)
        y = zero(ST)
        for j in eachindex(C.X)
            y += M.r₀[j] * C.X[j][k]
        end
        b[D*S+k] = sol.q[k] - y
    end
end


function update!(sol, params, int::GeometricIntegrator{<:CGVI}, DT)
    sol.q .= cache(int, DT).q̃
    sol.p .= cache(int, DT).p̃
end

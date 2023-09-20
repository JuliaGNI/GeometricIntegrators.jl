@doc raw"""
Continuous Galerkin Variational Integrator.

* `b`: weights of the quadrature rule
* `c`: nodes of the quadrature rule
* `x`: nodes of the basis
* `m`: mass matrix
* `a`: derivative matrix
* `r₀`: reconstruction coefficients at the beginning of the interval
* `r₁`: reconstruction coefficients at the end of the interval

"""
struct CGVI{T, NBASIS, NNODES, NDOF, basisType <: Basis{T}} <: LODEMethod
    basis::basisType
    quadrature::QuadratureRule{T,NNODES}

    b::SVector{NNODES,T}
    c::SVector{NNODES,T}

    x::SVector{NBASIS,T}

    m::SMatrix{NNODES, NBASIS, T, NDOF}
    a::SMatrix{NNODES, NBASIS, T, NDOF}

    r₀::SVector{NBASIS,T}
    r₁::SVector{NBASIS,T}

    function CGVI(basis::Basis{T}, quadrature::QuadratureRule{T}) where {T}
        # get number of quadrature nodes and number of basis functions
        NNODES = QuadratureRules.nnodes(quadrature)
        NBASIS = CompactBasisFunctions.nbasis(basis)

        # get quadrature nodes and weights
        quad_weights = QuadratureRules.weights(quadrature)
        quad_nodes = QuadratureRules.nodes(quadrature)

        # compute coefficients
        r₀ = zeros(T, NBASIS)
        r₁ = zeros(T, NBASIS)
        m  = zeros(T, NNODES, NBASIS)
        a  = zeros(T, NNODES, NBASIS)

        for i in eachindex(basis)
            r₀[i] = basis[zero(T), i]
            r₁[i] = basis[one(T), i]
            for j in eachindex(quad_nodes)
                m[j,i] = basis[quad_nodes[j], i]
                a[j,i] = basis'[quad_nodes[j], i]
            end
        end

        new{T, NBASIS, NNODES, NBASIS * NNODES, typeof(basis)}(basis, quadrature, quad_weights, quad_nodes, CompactBasisFunctions.grid(basis), m, a, r₀, r₁)
    end
end

basis(method::CGVI) = method.basis
quadrature(method::CGVI) = method.quadrature

nbasis(::CGVI{T,NB,NN}) where {T,NB,NN} = NB
nnodes(::CGVI{T,NB,NN}) where {T,NB,NN} = NN

isexplicit(::Union{CGVI, Type{<:CGVI}}) = false
isimplicit(::Union{CGVI, Type{<:CGVI}}) = true
issymmetric(::Union{CGVI, Type{<:CGVI}}) = missing
issymplectic(::Union{CGVI, Type{<:CGVI}}) = true

isiodemethod(::Union{CGVI, Type{<:CGVI}}) = true

default_solver(::CGVI) = Newton()
default_iguess(::CGVI) = HermiteExtrapolation()

function Base.show(io::IO, method::CGVI)
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


struct CGVICache{ST,D,S,R} <: IODEIntegratorCache{ST,D}
    x::Vector{ST}

    q̄::Vector{ST}
    p̄::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    s̃::Vector{ST}

    X::Vector{Vector{ST}}
    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}


    function CGVICache{ST,D,S,R}() where {ST,D,S,R}
        x = zeros(ST, D*(S+1))
        
        q̄ = zeros(ST,D)
        p̄ = zeros(ST,D)

        # create temporary vectors
        q̃ = zeros(ST,D)
        p̃ = zeros(ST,D)
        ṽ = zeros(ST,D)
        f̃ = zeros(ST,D)
        s̃ = zeros(ST,D)

        # create internal stage vectors
        X = create_internal_stage_vector(ST,D,S)
        Q = create_internal_stage_vector(ST,D,R)
        P = create_internal_stage_vector(ST,D,R)
        V = create_internal_stage_vector(ST,D,R)
        F = create_internal_stage_vector(ST,D,R)

        new(x, q̄, p̄, q̃, p̃, ṽ, f̃, s̃, X, Q, P, V, F)
    end
end

function reset!(cache::CGVICache, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

nlsolution(cache::CGVICache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::CGVI; kwargs...) where {ST}
    CGVICache{ST, ndims(problem), nbasis(method), nnodes(method)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::CGVI) = CGVICache{ST, ndims(problem), nbasis(method), nnodes(method)}


function initial_guess!(int::GeometricIntegrator{<:CGVI})
    # set some local variables for convenience
    local D = ndims(int)
    local S = nbasis(method(int))
    local x = nlsolution(int)

    for i in eachindex(basis(method(int)))
        initialguess!(method(int).x[i], cache(int).q̃, cache(int).p̃, solstep(int), problem(int), iguess(int))

        for k in 1:D
            x[D*(i-1)+k] = cache(int).q̃[k]
        end
    end

    initialguess!(one(timestep(int)), cache(int).q̃, cache(int).p̃, solstep(int), problem(int), iguess(int))

    for k in 1:D
        x[D*S+k] = cache(int).p̃[k]
    end
end


function components!(x::AbstractVector{ST}, int::GeometricIntegrator{<:CGVI}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local S = nbasis(method(int))
    local q = cache(int, ST).q̃
    local p = cache(int, ST).p̃
    local Q = cache(int, ST).Q
    local V = cache(int, ST).V
    local P = cache(int, ST).P
    local F = cache(int, ST).F
    local X = cache(int, ST).X


    # copy x to X
    for i in eachindex(X)
        for k in eachindex(X[i])
            X[i][k] = x[D*(i-1)+k]
        end
    end

    # copy x to p
    for k in eachindex(p)
        p[k] = x[D*S+k]
    end

    # compute Q
    for i in eachindex(Q)
        for k in eachindex(Q[i])
            y = zero(ST)
            for j in eachindex(X)
                y += method(int).m[i,j] * X[j][k]
            end
            Q[i][k] = y
        end
    end

    # compute q
    for k in eachindex(q)
        y = zero(ST)
        for i in eachindex(X)
            y += method(int).r₁[i] * X[i][k]
        end
        q[k] = y
    end

    # compute V
    for i in eachindex(V)
        for k in eachindex(V[i])
            y = zero(ST)
            for j in eachindex(X)
                y += method(int).a[i,j] * X[j][k]
            end
            V[i][k] = y / timestep(int)
        end
    end

    # compute P=ϑ(Q,V) and F=f(Q,V)
    for i in eachindex(Q,V,P,F)
        tᵢ = solstep(int).t + timestep(int) * method(int).c[i]
        equations(int).ϑ(P[i], tᵢ, Q[i], V[i], parameters(solstep(int)))
        equations(int).f(F[i], tᵢ, Q[i], V[i], parameters(solstep(int)))
    end
end


function residual!(b::Vector{ST}, int::GeometricIntegrator{<:CGVI}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local S = nbasis(method(int))
    local q̄ = cache(int, ST).q̄
    local p̄ = cache(int, ST).p̄
    local p̃ = cache(int, ST).p̃
    local P = cache(int, ST).P
    local F = cache(int, ST).F
    local X = cache(int, ST).X

    # compute b = - [(P-AF)]
    for i in eachindex(method(int).r₀, method(int).r₁)
        for k in eachindex(p̃, p̄)
            z = zero(ST)
            for j in eachindex(P,F)
                z += method(int).b[j] * method(int).m[j,i] * F[j][k] * timestep(int)
                z += method(int).b[j] * method(int).a[j,i] * P[j][k]
            end
            b[D*(i-1)+k] = (method(int).r₁[i] * p̃[k] - method(int).r₀[i] * p̄[k]) - z
        end
    end

    # compute b = - [(q-r₀Q)]
    for k in eachindex(q̄)
        y = zero(ST)
        for j in eachindex(X)
            y += method(int).r₀[j] * X[j][k]
        end
        b[D*S+k] = q̄[k] - y
    end
end


# Compute stages of Variational Partitioned Runge-Kutta methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:CGVI}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:CGVI}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int, DT).q̃
    solstep(int).p .= cache(int, DT).p̃
end


function integrate_step!(int::GeometricIntegrator{<:CGVI, <:AbstractProblemIODE})
    # copy previous solution from solstep to cache
    reset!(cache(int), current(solstep(int))...)

    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    # check_solver_status(int.solver.status, int.solver.params)

    # compute final update
    update!(nlsolution(int), int)
end

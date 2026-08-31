@doc raw"""
Variational partitioned Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of q
* `p̃`: initial guess of p
* `ṽ`: initial guess of v
* `f̃`: initial guess of f
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of q
* `P`: internal stages of p
* `V`: internal stages of v
* `F`: internal stages of f
* `Y`: integral of vector field of internal stages of q
* `Z`: integral of vector field of internal stages of p
"""
mutable struct VPRKCache{ST, S} <: IODEIntegratorCache{ST}
    x::Vector{ST}
    x̄::Vector{ST}

    λ::Vector{ST}
    λ̄::Vector{ST}

    q₋::Vector{ST}
    q̄₊::Vector{ST}
    p₋::Vector{ST}
    p̄₊::Vector{ST}

    u::Vector{ST}
    g::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    θ̃::Vector{ST}
    s̃::Vector{ST}

    μ::Vector{ST}
    ϕ::Vector{ST}
    v::Vector{ST}
    f::Vector{ST}
    y::Vector{ST}
    z::Vector{ST}

    Q::Vector{Vector{ST}}
    P::Vector{Vector{ST}}
    V::Vector{Vector{ST}}
    F::Vector{Vector{ST}}
    Λ::Vector{Vector{ST}}
    Φ::Vector{Vector{ST}}

    Y::Vector{Vector{ST}}
    Z::Vector{Vector{ST}}
    U::Vector{Vector{ST}}
    G::Vector{Vector{ST}}
    R::Vector{Vector{ST}}

    function VPRKCache{ST, S}(ics, n, projection::Bool = false, m = 0) where {ST, S}
        D = length(vec(ics.q))

        # create solver vector
        x = zeros(ST, n)
        x̄ = zeros(ST, m)

        # create temporary vectors
        q̃ = zeros(ST, D)
        p̃ = zeros(ST, D)
        ṽ = zeros(ST, D)
        f̃ = zeros(ST, D)
        θ̃ = zeros(ST, D)
        s̃ = zeros(ST, D)

        ϕ = zeros(ST, D)
        μ = zeros(ST, D)

        # create update vectors
        v = zeros(ST, D)
        f = zeros(ST, D)
        y = zeros(ST, D)
        z = zeros(ST, D)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        # projection vectors
        if projection
            λ = zeros(ST, D)
            λ̄ = zeros(ST, D)

            q₋ = zeros(ST, D)
            q̄₊ = zeros(ST, D)
            p₋ = zeros(ST, D)
            p̄₊ = zeros(ST, D)

            u = zeros(ST, D)
            g = zeros(ST, D)

            Λ = create_internal_stage_vector(ST, D, S)
            Φ = create_internal_stage_vector(ST, D, S)

            U = create_internal_stage_vector(ST, D, 2)
            G = create_internal_stage_vector(ST, D, 2)
            R = create_internal_stage_vector(ST, D, S)
        else
            λ = Vector{ST}()
            λ̄ = Vector{ST}()

            q₋ = Vector{ST}()
            q̄₊ = Vector{ST}()
            p₋ = Vector{ST}()
            p̄₊ = Vector{ST}()

            u = Vector{ST}()
            g = Vector{ST}()

            Λ = create_internal_stage_vector(ST, 0, 0)
            Φ = create_internal_stage_vector(ST, 0, 0)

            U = create_internal_stage_vector(ST, 0, 0)
            G = create_internal_stage_vector(ST, 0, 0)
            R = create_internal_stage_vector(ST, 0, 0)
        end

        new(x, x̄, λ, λ̄, q₋, q̄₊, p₋, p̄₊,
            u, g, q̃, p̃, ṽ, f̃, θ̃, s̃, ϕ, μ, v, f, y, z,
            Q, P, V, F, Λ, Φ, Y, Z, U, G, R)
    end
end

function VPRKCache(ics, ST, S, N, ::VPRKMethod)
    VPRKCache{ST, S}(ics, N, false)
end

# function IntegratorCacheVPRK(ics, ST, S, N, ::ProjectedVPRK)
#     IntegratorCacheVPRK{ST,S}(ics, N, true)
# end

nlsolution(cache::VPRKCache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::VPRKMethod; kwargs...) where {ST}
    S = nstages(method)
    N = solversize(method, problem)

    # if hasnullvector(method)
    #     N += ndims(problem)
    # end

    VPRKCache(initial_conditions(problem), ST, S, N, method; kwargs...)
end

@inline CacheType(ST, ::AbstractProblemIODE, method::VPRKMethod) = VPRKCache{
    ST, nstages(method)}

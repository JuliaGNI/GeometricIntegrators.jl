"""
Variational partitioned Runge-Kutta integrator cache.

### Fields

* `n`: time step number
* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `p`: current solution of p
* `p̅`: previous solution of p
* `v`: vector field of q
* `v̅`: vector field of q̅
* `f`: vector field of p
* `f̅`: vector field of p̅
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
mutable struct IntegratorCacheVPRK{ST,TT,D,S} <: IODEIntegratorCache{ST,D}
    n::Int
    t::TT
    t̅::TT

    q::Vector{ST}
    q̅::Vector{ST}
    p::Vector{ST}
    p̅::Vector{ST}
    λ::Vector{ST}
    λ̅::Vector{ST}

    qₑᵣᵣ::Vector{ST}
    pₑᵣᵣ::Vector{ST}

    q₋::Vector{ST}
    q̅₊::Vector{ST}
    p₋::Vector{ST}
    p̅₊::Vector{ST}

    v::Vector{ST}
    v̅::Vector{ST}
    f::Vector{ST}
    f̅::Vector{ST}

    u::Vector{ST}
    g::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    s̃::Vector{ST}

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

    function IntegratorCacheVPRK{ST,TT,D,S}(projection::Bool=false) where {ST,TT,D,S}
        # create solution vectors
        q = zeros(ST,D)
        q̅ = zeros(ST,D)
        p = zeros(ST,D)
        p̅ = zeros(ST,D)

        # create error vectors
        qₑᵣᵣ = zeros(ST,D)
        pₑᵣᵣ = zeros(ST,D)

        # create update vectors
        v = zeros(ST,D)
        v̅ = zeros(ST,D)
        f = zeros(ST,D)
        f̅ = zeros(ST,D)

        # create temporary vectors
        q̃ = zeros(ST,D)
        p̃ = zeros(ST,D)
        ṽ = zeros(ST,D)
        f̃ = zeros(ST,D)
        s̃ = zeros(ST,D)

        # create internal stage vectors
        Q = create_internal_stage_vector(ST, D, S)
        P = create_internal_stage_vector(ST, D, S)
        V = create_internal_stage_vector(ST, D, S)
        F = create_internal_stage_vector(ST, D, S)
        Y = create_internal_stage_vector(ST, D, S)
        Z = create_internal_stage_vector(ST, D, S)

        # projection vectors
        if projection
            λ = zeros(ST,D)
            λ̅ = zeros(ST,D)

            q₋= zeros(ST,D)
            q̅₊= zeros(ST,D)
            p₋= zeros(ST,D)
            p̅₊= zeros(ST,D)

            u = zeros(ST,D)
            g = zeros(ST,D)

            Λ = create_internal_stage_vector(ST, D, S)
            Φ = create_internal_stage_vector(ST, D, S)

            U = create_internal_stage_vector(ST, D, 2)
            G = create_internal_stage_vector(ST, D, 2)
            R = create_internal_stage_vector(ST, D, S)
        else
            λ = Vector{ST}()
            λ̅ = Vector{ST}()

            q₋= Vector{ST}()
            q̅₊= Vector{ST}()
            p₋= Vector{ST}()
            p̅₊= Vector{ST}()

            u = Vector{ST}()
            g = Vector{ST}()

            Λ = create_internal_stage_vector(ST, 0, 0)
            Φ = create_internal_stage_vector(ST, 0, 0)

            U = create_internal_stage_vector(ST, 0, 0)
            G = create_internal_stage_vector(ST, 0, 0)
            R = create_internal_stage_vector(ST, 0, 0)
        end

        new(0, zero(TT), zero(TT),
            q, q̅, p, p̅, λ, λ̅,
            qₑᵣᵣ, pₑᵣᵣ,
            q₋, q̅₊, p₋, p̅₊,
            v, v̅, f, f̅, u, g, q̃, p̃, ṽ, f̃, s̃,
            Q, P, V, F, Λ, Φ, Y, Z, U, G, R)
    end
end

function IntegratorCacheVPRK(ST,TT,D,S)
    IntegratorCacheVPRK{ST,TT,D,S}(false)
end

function IntegratorCacheVPRKwProjection(ST,TT,D,S)
    IntegratorCacheVPRK{ST,TT,D,S}(true)
end

function update_params!(params::AbstractParametersVPRK, cache::IntegratorCacheVPRK)
    # set time for nonlinear solver and copy previous solution
    params.t̅  = cache.t
    params.q̅ .= cache.q
    params.p̅ .= cache.p
end

function update_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::IntegratorCacheVPRK{DT,TT}) where {DT,TT}
    update_solution!(cache.q, cache.qₑᵣᵣ, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(cache.p, cache.pₑᵣᵣ, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))
end

function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::IntegratorCacheVPRK, R::Vector{TT}) where {DT,TT}
    update_solution!(cache.q, cache.qₑᵣᵣ, cache.U, R, timestep(int))
    update_solution!(cache.p, cache.pₑᵣᵣ, cache.G, R, timestep(int))
end

function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, cache::IntegratorCacheVPRK, RU::Vector{TT}, RG::Vector{TT}) where {DT,TT}
    update_solution!(cache.q, cache.qₑᵣᵣ, cache.U, RU, timestep(int))
    update_solution!(cache.p, cache.pₑᵣᵣ, cache.G, RG, timestep(int))
end

function CommonFunctions.set_solution!(cache::IntegratorCacheVPRK, sol, n=0)
    t, q, p = sol
    cache.n  = n
    cache.t  = t
    cache.q .= q
    cache.p .= p
    cache.v .= 0
    cache.f .= 0
end


function create_integrator_cache(int::AbstractIntegratorVPRK{DT,TT}) where {DT,TT}
    IntegratorCacheVPRK(DT, TT, ndims(int), nstages(int))
end


function create_integrator_cache(int::AbstractIntegratorVPRKwProjection{DT,TT}) where {DT,TT}
    IntegratorCacheVPRKwProjection(DT, TT, ndims(int), nstages(int))
end


function initialize!(int::AbstractIntegratorVPRK{DT}, cache::IntegratorCacheVPRK) where {DT,TT}
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end

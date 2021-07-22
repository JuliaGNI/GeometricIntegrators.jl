@doc raw"""
Cache of a Specialised Partitioned Additive Runge-Kutta integrator.

### Fields

* `n`: time step number
* `t`: time of current time step
* `t̄`: time of previous time step
* `q`: current solution of q
* `q̄`: previous solution of q
* `p`: current solution of p
* `p̄`: previous solution of p
* `v`: vector field of q
* `v̄`: vector field of q̄
* `f`: vector field of p
* `f̄`: vector field of p̄
* `q̃`: initial guess of q
* `p̃`: initial guess of p
* `ṽ`: initial guess of v
* `f̃`: initial guess of f
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of q
* `P`: internal stages of p
* `V`: internal stages of v
* `F`: internal stages of f
* `Y`: vector field of internal stages of q
* `Z`: vector field of internal stages of p
"""
mutable struct IntegratorCacheSPARK{ST,D,S,R} <: IDAEIntegratorCache{ST,D}
    q::Vector{ST}
    q̄::Vector{ST}
    p::Vector{ST}
    p̄::Vector{ST}
    λ::Vector{ST}
    λ̄::Vector{ST}
    μ::Vector{ST}
    μ̅::Vector{ST}

    v::Vector{ST}
    v̄::Vector{ST}
    f::Vector{ST}
    f̄::Vector{ST}
    u::Vector{ST}
    ū::Vector{ST}
    g::Vector{ST}
    ḡ::Vector{ST}
    ϕ::Vector{ST}
    ϕ̅::Vector{ST}

    q̃::Vector{ST}
    p̃::Vector{ST}
    ṽ::Vector{ST}
    f̃::Vector{ST}
    ϕ̃::Vector{ST}
    s̃::Vector{ST}

    Qi::Vector{Vector{ST}}
    Pi::Vector{Vector{ST}}
    Vi::Vector{Vector{ST}}
    Fi::Vector{Vector{ST}}
    Gi::Vector{Vector{ST}}
    Hi::Vector{Vector{ST}}
    Yi::Vector{Vector{ST}}
    Zi::Vector{Vector{ST}}
    Φi::Vector{Vector{ST}}
    Ψi::Vector{Vector{ST}}

    Qp::Vector{Vector{ST}}
    Pp::Vector{Vector{ST}}
    Vp::Vector{Vector{ST}}
    Λp::Vector{Vector{ST}}
    Up::Vector{Vector{ST}}
    Fp::Vector{Vector{ST}}
    Gp::Vector{Vector{ST}}
    G̅p::Vector{Vector{ST}}
    Hp::Vector{Vector{ST}}
    Yp::Vector{Vector{ST}}
    Zp::Vector{Vector{ST}}
    Φp::Vector{Vector{ST}}
    Ψp::Vector{Vector{ST}}

    function IntegratorCacheSPARK{ST,D,S,R}() where {ST,D,S,R}
        q = zeros(ST,D)
        q̄ = zeros(ST,D)
        p = zeros(ST,D)
        p̄ = zeros(ST,D)
        λ = zeros(ST,D)
        λ̄= zeros(ST,D)
        μ = zeros(ST,D)
        μ̅ = zeros(ST,D)

        # create update vectors
        v = zeros(ST,D)
        v̄ = zeros(ST,D)
        f = zeros(ST,D)
        f̄ = zeros(ST,D)
        u = zeros(ST,D)
        ū = zeros(ST,D)
        g = zeros(ST,D)
        ḡ = zeros(ST,D)
        ϕ = zeros(ST,D)
        ϕ̅ = zeros(ST,D)

        # create temporary vectors
        q̃ = zeros(ST,D)
        p̃ = zeros(ST,D)
        ṽ = zeros(ST,D)
        f̃ = zeros(ST,D)
        ϕ̃ = zeros(ST,D)
        s̃ = zeros(ST,D)

        # create internal stage vectors
        Qi = create_internal_stage_vector(ST, D, S)
        Pi = create_internal_stage_vector(ST, D, S)
        Vi = create_internal_stage_vector(ST, D, S)
        Fi = create_internal_stage_vector(ST, D, S)
        Gi = create_internal_stage_vector(ST, D, S)
        Hi = create_internal_stage_vector(ST, D, S)
        Yi = create_internal_stage_vector(ST, D, S)
        Zi = create_internal_stage_vector(ST, D, S)
        Φi = create_internal_stage_vector(ST, D, S)
        Ψi = create_internal_stage_vector(ST, D, S)

        Qp = create_internal_stage_vector(ST, D, R)
        Pp = create_internal_stage_vector(ST, D, R)
        Vp = create_internal_stage_vector(ST, D, R)
        Λp = create_internal_stage_vector(ST, D, R)
        Up = create_internal_stage_vector(ST, D, R)
        Fp = create_internal_stage_vector(ST, D, R)
        Gp = create_internal_stage_vector(ST, D, R)
        G̅p = create_internal_stage_vector(ST, D, R)
        Hp = create_internal_stage_vector(ST, D, R)
        Yp = create_internal_stage_vector(ST, D, R)
        Zp = create_internal_stage_vector(ST, D, R)
        Φp = create_internal_stage_vector(ST, D, R)
        Ψp = create_internal_stage_vector(ST, D, R)

        new(q, q̄, p, p̄, λ, λ̄, μ, μ̅,
            v, v̄, f, f̄, u, ū, g, ḡ, ϕ, ϕ̅,
            q̃, p̃, ṽ, f̃, ϕ̃, s̃,
            Qi, Pi, Vi, Fi, Gi, Hi, Yi, Zi, Φi, Ψi,
            Qp, Pp, Vp, Λp, Up, Fp, Gp, G̅p, Hp, Yp, Zp, Φp, Ψp)
    end
end

function IntegratorCache(int::AbstractIntegratorSPARK{DT,TT}) where {DT,TT}
    IntegratorCacheSPARK{DT, ndims(int), nstages(int), pstages(int)}()
end

function GeometricBase.reset!(cache::IntegratorCacheSPARK, Δt)
    cache.t̄  = cache.t
    cache.q̄ .= cache.q
    cache.p̄ .= cache.p
    cache.λ̄.= cache.λ
    cache.v̄ .= cache.v
    cache.f̄ .= cache.f
    cache.ū .= cache.u
    cache.ḡ .= cache.g
    cache.ϕ̅ .= cache.ϕ
end

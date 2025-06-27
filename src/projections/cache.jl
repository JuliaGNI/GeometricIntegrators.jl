
mutable struct ProjectionCache{DT,TT,PT,D,M,N} <: IODEIntegratorCache{DT,D}
    t::TT

    x::Vector{DT}
    x̄::SubArray{DT, 1, Vector{DT}, Tuple{UnitRange{Int}}, true}
    x̃::SubArray{DT, 1, Vector{DT}, Tuple{UnitRange{Int}}, true}

    q::Vector{DT}
    p::Vector{DT}
    v::Vector{DT}
    f::Vector{DT}

    q̄::Vector{DT}
    q̃::Vector{DT}
    p̃::Vector{DT}
    ṽ::Vector{DT}
    f̃::Vector{DT}

    λ::Vector{DT}
    ϑ::Vector{DT}
    ϕ::Vector{DT}
    u::Vector{DT}
    g::Vector{DT}

    U::Vector{Vector{DT}}
    G::Vector{Vector{DT}}

    function ProjectionCache{DT}(problem::EquationProblem, method::ProjectedMethod) where {DT}
        D = ndims(problem)
        M = nconstraints(problem)
        N = solversize(problem, parent(method))

        TT = timetype(problem)
        t = initialtime(problem)

        x = zeros(DT, N+D+M)
        x̄ = @view x[1:N]
        x̃ = @view x[N+1:N+D+M]

        q = zeros(DT, D)
        p = zeros(DT, D)
        v = zeros(DT, D)
        f = zeros(DT, D)

        q̄ = zeros(DT, D)
        q̃ = zeros(DT, D)
        p̃ = zeros(DT, D)
        ṽ = zeros(DT, D)
        f̃ = zeros(DT, D)

        λ = zeros(DT, M)
        ϑ = zeros(DT, M)
        ϕ = zeros(DT, M)
        u = zeros(DT, D)
        g = zeros(DT, D)

        U = [zeros(DT, D), zeros(DT, D)]
        G = [zeros(DT, D), zeros(DT, D)]

        new{DT, TT, typeof(problem), D, M, N}(t, x, x̄, x̃, q, p, v, f, q̄, q̃, p̃, ṽ, f̃, λ, ϑ, ϕ, u, g, U, G)
    end
end

ProjectionCache(problem::EquationProblem, method::ProjectedMethod) = ProjectionCache{datatype(problem)}(problem, method)

Base.ndims(::ProjectionCache{DT,TT,PT,D,M,N}) where {DT,TT,PT,D,M,N} = D
nconstraints(::ProjectionCache{DT,TT,PT,D,M,N}) where {DT,TT,PT,D,M,N} = M
solversize(::ProjectionCache{DT,TT,PT,D,M,N}) where {DT,TT,PT,D,M,N} = N

nlsolution(cache::ProjectionCache) = cache.x

function split_nlsolution(cache::ProjectionCache{DT,TT,PT,D,M,N}) where {DT,TT,PT,D,M,N}
    x = nlsolution(cache)
    x̄ = @view x[1:N]
    x̃ = @view x[N+1:N+D+M]

    return (x̄, x̃)
end


function current(cache::ProjectionCache{DT, TT, <:DAEProblem}) where {DT, TT}
    (t = cache.t, q = cache.q)
end

function current(cache::ProjectionCache{DT, TT, <:Union{PDAEProblem, HDAEProblem}}) where {DT, TT}
    (t = cache.t, q = cache.q, p = cache.p)
end

function current(cache::ProjectionCache{DT, TT, <:Union{IODEProblem, LODEProblem}}) where {DT, TT}
    (t = cache.t, q = cache.q, v = cache.v, p = cache.p)
end

function reset!(cache::ProjectionCache{DT, TT, <:DAEProblem}, sol) where {DT, TT}
    cache.t = sol.t
    copyto!(cache.q, sol.q)
    return cache
end

function reset!(cache::ProjectionCache{DT, TT, <:Union{PDAEProblem, HDAEProblem}}, sol) where {DT, TT}
    cache.t = sol.t
    copyto!(cache.q, sol.q)
    copyto!(cache.p, sol.p)
    return cache
end

function reset!(cache::ProjectionCache{DT, TT, <:Union{IODEProblem, LODEProblem}}, sol) where {DT, TT}
    cache.t = sol.t
    copyto!(cache.q, sol.q)
    copyto!(cache.v, sol.v)
    copyto!(cache.p, sol.p)
    return cache
end

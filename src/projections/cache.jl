
struct ProjectionCache{DT,D,M,N} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}
    x̄::SubArray{DT, 1, Vector{DT}, Tuple{UnitRange{Int}}, true}
    x̃::SubArray{DT, 1, Vector{DT}, Tuple{UnitRange{Int}}, true}

    q::Vector{DT}
    q̃::Vector{DT}

    p::Vector{DT}

    Δq::Vector{DT}
    Δp::Vector{DT}

    v::Vector{DT}
    f::Vector{DT}

    λ::Vector{DT}
    ϕ::Vector{DT}
    u::Vector{DT}
    g::Vector{DT}

    U::Vector{Vector{DT}}
    G::Vector{Vector{DT}}

    function ProjectionCache{DT,D,M,N}() where {DT,D,M,N}
        x = zeros(DT, N+D+M)
        x̄ = @view x[1:N]
        x̃ = @view x[N+1:N+D+M]

        q = zeros(DT, D)
        q̃ = zeros(DT, D)

        p = zeros(DT, D)

        Δq = zeros(DT, D)
        Δp = zeros(DT, D)

        v = zeros(DT, D)
        f = zeros(DT, D)
        λ = zeros(DT, M)
        ϕ = zeros(DT, M)
        u = zeros(DT, D)
        g = zeros(DT, D)
        U = [zeros(DT, D), zeros(DT, D)]
        G = [zeros(DT, D), zeros(DT, D)]

        new(x, x̄, x̃, q, q̃, p, Δq, Δp, v, f, λ, ϕ, u, g, U, G)
    end
end

Base.ndims(::ProjectionCache{DT,D,M,N}) where {DT,D,M,N} = D
nconstraints(::ProjectionCache{DT,D,M,N}) where {DT,D,M,N} = M
solversize(::ProjectionCache{DT,D,M,N}) where {DT,D,M,N} = N

nlsolution(cache::ProjectionCache) = cache.x

function split_nlsolution(cache::ProjectionCache{DT,D,M,N}) where {DT,D,M,N}
    x = nlsolution(cache)
    x̄ = @view x[1:N]
    x̃ = @view x[N+1:N+D+M]

    return (x̄, x̃)
end


function current(cache::ProjectionCache, solstep::Union{SolutionStepODE,SolutionStepDAE})
    (t = solstep.t, q = cache.q)
end

function current(cache::ProjectionCache, solstep::Union{SolutionStepPODE,SolutionStepPDAE})
    (t = solstep.t, q = cache.q, p = cache.p)
end

function reset!(cache::ProjectionCache, t, q)
    copyto!(cache.q, q)
end

function reset!(cache::ProjectionCache, t, q, p, λ = missing)
    copyto!(cache.q, q)
    copyto!(cache.p, p)
end

function update!(cache::ProjectionCache, V::AbstractVector, tableau::Tableau, Δt::Number)
    update_vector!(cache.Δq, V, tableau, Δt)
    cache.q .+= cache.Δq
end

function update!(cache::ProjectionCache, V::AbstractVector, F::AbstractVector, tableau::Tableau, Δt::Number)
    update_vector!(cache.Δq, V, tableau, Δt)
    update_vector!(cache.Δp, F, tableau, Δt)
    cache.q .+= cache.Δq
    cache.p .+= cache.Δp
end

function update!(cache::ProjectionCache, V::AbstractVector, F::AbstractVector, tableau::PartitionedTableau, Δt::Number)
    update_vector!(cache.Δq, V, tableau.q, Δt)
    update_vector!(cache.Δp, F, tableau.p, Δt)
    cache.q .+= cache.Δq
    cache.p .+= cache.Δp
end

function update!(cache::ProjectionCache, subcache::Union{ODEIntegratorCache,DAEIntegratorCache}, tableau::Union{Tableau,PartitionedTableau}, Δt::Number)
    update!(cache, subcache.V, tableau, Δt)
end

function update!(cache::ProjectionCache, subcache::Union{IODEIntegratorCache,PODEIntegratorCache}, tableau::Union{Tableau,PartitionedTableau}, Δt::Number)
    update!(cache, subcache.V, subcache.F, tableau, Δt)
end

function project!(cache::ProjectionCache, U::AbstractVector, G::AbstractVector, Δt::Number)
    cache.q .+= Δt .* U
    cache.p .+= Δt .* G
end

function project!(solstep::SolutionStep, problem::EquationProblem, method::ProjectionMethod, cache::ProjectionCache)
    project(solstep, problem, method, cache.U, cache.G, cache.λ)
end


"Trapezoidal Degenerate Variational Integrator."
struct CTDVI <: DVIMethod end

order(::Union{CTDVI, Type{CTDVI}}) = 2

issymmetric(::Union{CTDVI, Type{<:CTDVI}}) = true


@doc raw"""
Degenerate variational integrator cache.

### Fields

* `q`: internal stages of solution
* `v`: internal stages of vector field
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
"""
struct CTDVICache{DT,D} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}

    q̄::Vector{DT}
    p̄::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    q⁻::Vector{DT}
    θ⁻::Vector{DT}
    f⁻::Vector{DT}

    q⁺::Vector{DT}
    θ⁺::Vector{DT}
    f⁺::Vector{DT}

    function CTDVICache{DT,D}() where {DT,D}
        new(zeros(DT,2D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function reset!(cache::CTDVICache, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

nlsolution(cache::CTDVICache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::CTDVI; kwargs...) where {ST}
    CTDVICache{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, ::CTDVI) = CTDVICache{ST, ndims(problem)}


function Base.show(io::IO, int::GeometricIntegrator{<:CTDVI})
    print(io, "\nTrapezoidal Degenerate Variational Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function initial_guess!(int::GeometricIntegrator{<:CTDVI})
    # set some local variables for convenience
    local D = ndims(int)
    local x = nlsolution(int)

    # compute initial guess for solution q(n+1)
    initialguess!(solstep(int).t, cache(int).q, cache(int).θ, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

    x[1:D] .= cache(int).q

    # compute initial guess for solution q(n+1/2)
    initialguess!((solstep(int).t + solstep(int).t̄)/2, cache(int).q, cache(int).θ, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

    offset_v = D
    offset_x = D + div(D,2)
    for k in 1:div(D,2)
        cache(int).x[offset_v+k] = cache(int).v[k]              # v¹(n+1/2)
        cache(int).x[offset_x+k] = cache(int).q[div(D,2)+k]     # q²(n+1/2)
    end
end


function components!(x::Vector{ST}, int::GeometricIntegrator{<:CTDVI}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local t⁻ = solstep(int).t̄
    local t⁺ = solstep(int).t
    
    # copy x to q
    cache(int,ST).q .= x[1:D]

    # copy x to q⁻, q⁺ and v
    for k in 1:div(D,2)
        cache(int,ST).q⁻[k] = cache(int).q̄[k]
        cache(int,ST).q⁺[k] = cache(int,ST).q[k]

        cache(int,ST).q⁻[div(D,2)+k] = x[D+div(D,2)+k]
        cache(int,ST).q⁺[div(D,2)+k] = x[D+div(D,2)+k]

        cache(int,ST).v[k]          = x[D+k]
        cache(int,ST).v[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    equations(int).f(cache(int,ST).f⁻, t⁻, cache(int,ST).q⁻, cache(int,ST).v, parameters(solstep(int)))
    equations(int).f(cache(int,ST).f⁺, t⁺, cache(int,ST).q⁺, cache(int,ST).v, parameters(solstep(int)))
 
    # compute Θ = ϑ(q,v)
    equations(int).ϑ(cache(int,ST).θ⁻, t⁻, cache(int,ST).q⁻, cache(int,ST).v, parameters(solstep(int)))
    equations(int).ϑ(cache(int,ST).θ⁺, t⁺, cache(int,ST).q⁺, cache(int,ST).v, parameters(solstep(int)))
    equations(int).ϑ(cache(int,ST).θ,  t⁺, cache(int,ST).q,  cache(int,ST).v, parameters(solstep(int)))
end


function residual!(b::Vector{ST}, x::Vector{ST}, int::GeometricIntegrator{<:CTDVI}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute b
    b[1:D] .= (cache(int,ST).θ⁻ .+ cache(int,ST).θ⁺) ./ 2 .- cache(int).p̄ .- timestep(int) .* cache(int,ST).f⁻ ./ 2
    
    for k in 1:div(D,2)
        b[D+k]          = cache(int,ST).q[k] - cache(int).q̄[k] - timestep(int) * cache(int,ST).v[k]
        b[D+div(D,2)+k] = cache(int,ST).θ[k] - cache(int).p̄[k] - timestep(int) * (cache(int,ST).f⁻[k] + cache(int,ST).f⁺[k]) / 2
    end
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:CTDVI}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int, DT).q
    solstep(int).p .= cache(int, DT).θ
end

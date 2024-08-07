
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


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CTDVI})
    # set some local variables for convenience
    local D = ndims(int)
    local x = nlsolution(int)

    # compute initial guess for solution q(n+1)
    soltmp = (
        t = sol.t,
        q = cache(int).q,
        p = cache(int).θ,
        v = cache(int).v,
        f = cache(int).f,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    x[1:D] .= cache(int).q

    # compute initial guess for solution q(n+1/2)
    soltmp = (
        t = (sol.t + history.t[1]) / 2,
        q = cache(int).q,
        p = cache(int).θ,
        v = cache(int).v,
        f = cache(int).f,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    offset_v = D
    offset_x = D + div(D,2)
    for k in 1:div(D,2)
        cache(int).x[offset_v+k] = cache(int).v[k]              # v¹(n+1/2)
        cache(int).x[offset_x+k] = cache(int).q[div(D,2)+k]     # q²(n+1/2)
    end
end


function components!(x::Vector{ST}, sol, params, int::GeometricIntegrator{<:CTDVI}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local t⁻ = sol.t - timestep(int)
    local t⁺ = sol.t
    
    # copy x to q
    cache(int,ST).q .= x[1:D]

    # copy x to q⁻, q⁺ and v
    for k in 1:div(D,2)
        cache(int,ST).q⁻[k] = sol.q[k]
        cache(int,ST).q⁺[k] = cache(int,ST).q[k]

        cache(int,ST).q⁻[div(D,2)+k] = x[D+div(D,2)+k]
        cache(int,ST).q⁺[div(D,2)+k] = x[D+div(D,2)+k]

        cache(int,ST).v[k]          = x[D+k]
        cache(int,ST).v[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    equations(int).f(cache(int,ST).f⁻, t⁻, cache(int,ST).q⁻, cache(int,ST).v, params)
    equations(int).f(cache(int,ST).f⁺, t⁺, cache(int,ST).q⁺, cache(int,ST).v, params)
 
    # compute Θ = ϑ(q,v)
    equations(int).ϑ(cache(int,ST).θ⁻, t⁻, cache(int,ST).q⁻, cache(int,ST).v, params)
    equations(int).ϑ(cache(int,ST).θ⁺, t⁺, cache(int,ST).q⁺, cache(int,ST).v, params)
    equations(int).ϑ(cache(int,ST).θ,  t⁺, cache(int,ST).q,  cache(int,ST).v, params)
end


function residual!(b::Vector{ST}, sol, params, int::GeometricIntegrator{<:CTDVI}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # compute b
    b[1:D] .= (cache(int,ST).θ⁻ .+ cache(int,ST).θ⁺) ./ 2 .- sol.p .- timestep(int) .* cache(int,ST).f⁻ ./ 2
    
    for k in 1:div(D,2)
        b[D+k]          = cache(int,ST).q[k] - sol.q[k] - timestep(int) * cache(int,ST).v[k]
        b[D+div(D,2)+k] = cache(int,ST).θ[k] - sol.p[k] - timestep(int) * (cache(int,ST).f⁻[k] + cache(int,ST).f⁺[k]) / 2
    end
end


function update!(sol, params, int::GeometricIntegrator{<:CTDVI}, DT)
    # compute final update
    sol.q .= cache(int, DT).q
    sol.p .= cache(int, DT).θ
end

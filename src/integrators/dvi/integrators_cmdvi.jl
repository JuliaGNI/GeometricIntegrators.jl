
"Midpoint Degenerate Variational Integrator."
struct CMDVI <: DVIMethod end

order(::Union{CMDVI, Type{CMDVI}}) = 2
issymmetric(::Union{CMDVI, Type{<:CMDVI}}) = true


function Base.show(io::IO, int::GeometricIntegrator{<:CMDVI})
    print(io, "\nMidpoint Degenerate Variational Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function initial_guess!(int::GeometricIntegrator{<:CMDVI})
    # set some local variables for convenience
    local D = ndims(int)
    local x = nlsolution(int)

    # compute initial guess for solution q(n+1)
    initialguess!(solstep(int).t, cache(int).q, cache(int).p, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

    x[1:D] .= cache(int).q

    # compute initial guess for solution q(n+1/2)
    initialguess!((solstep(int).t + solstep(int).t̄)/2, cache(int).q, cache(int).p, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

    offset_v = D
    offset_x = D + div(D,2)
    for k in 1:div(D,2)
        x[offset_v+k] = cache(int).v[k]              # v¹(n+1/2)
        x[offset_x+k] = cache(int).q[div(D,2)+k]     # q²(n+1/2)
    end
end


function components!(x::Vector{ST}, int::GeometricIntegrator{<:CMDVI}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local t̃ = solstep(int).t̄ + timestep(int) / 2
    local t = solstep(int).t
    
    # copy x to q
    cache(int,ST).q .= x[1:D]

    # copy x to q⁻, q⁺ and v
    for k in 1:div(D,2)
        cache(int,ST).q̃[k] = (cache(int).q̄[k] + cache(int,ST).q[k]) / 2
        cache(int,ST).q̃[div(D,2)+k] = x[D+div(D,2)+k]

        cache(int,ST).v[k]          = x[D+k]
        cache(int,ST).v[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    equations(int).f(cache(int,ST).f̃, t̃, cache(int,ST).q̃, cache(int,ST).v)
 
    # compute Θ = ϑ(q,v)
    equations(int).ϑ(cache(int,ST).p̃, t̃, cache(int,ST).q̃, cache(int,ST).v)
    equations(int).ϑ(cache(int,ST).p, t, cache(int,ST).q, cache(int,ST).v)
end


function residual!(b::Vector{ST}, x::Vector{ST}, int::GeometricIntegrator{<:CMDVI}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute b
    b[1:D] .= cache(int,ST).p̃ .- cache(int).p̄ .- timestep(int) .* cache(int,ST).f̃ ./ 2
    
    for k in 1:div(D,2)
        b[D+k]          = cache(int,ST).q[k] - cache(int).q̄[k] - timestep(int) * cache(int,ST).v[k]
        b[D+div(D,2)+k] = cache(int,ST).p[k] - cache(int).p̄[k] - timestep(int) * cache(int,ST).f̃[k]
    end
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:CMDVI}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int, DT).q
    solstep(int).p .= cache(int, DT).p
end


"Midpoint Degenerate Variational Integrator."
struct CMDVI <: DVIMethod end

order(::Union{CMDVI, Type{CMDVI}}) = 2
issymmetric(::Union{CMDVI, Type{<:CMDVI}}) = true


function Base.show(io::IO, int::GeometricIntegrator{<:CMDVI})
    print(io, "\nMidpoint Degenerate Variational Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CMDVI})
    # set some local variables for convenience
    local D = ndims(int)
    local x = nlsolution(int)

    # compute initial guess for solution q(n+1)
    soltmp = (
        t = sol.t,
        q = cache(int).q,
        p = cache(int).p,
        v = cache(int).v,
        f = cache(int).f,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    x[1:D] .= cache(int).q

    # compute initial guess for solution q(n+1/2)
    soltmp = (
        t = (sol.t + history.t[1]) / 2,
        q = cache(int).q,
        p = cache(int).p,
        v = cache(int).v,
        f = cache(int).f,
    )
    solutionstep!(soltmp, history, problem(int), iguess(int))

    offset_v = D
    offset_x = D + div(D,2)
    for k in 1:div(D,2)
        x[offset_v+k] = cache(int).v[k]              # v¹(n+1/2)
        x[offset_x+k] = cache(int).q[div(D,2)+k]     # q²(n+1/2)
    end
end


function components!(x::Vector{ST}, sol, params, int::GeometricIntegrator{<:CMDVI}) where {ST}
    # set some local variables for convenience and clarity
    local D = ndims(int)
    local t̃ = sol.t - timestep(int) / 2
    local t = sol.t
    
    # copy x to q
    cache(int,ST).q .= x[1:D]

    # copy x to q⁻, q⁺ and v
    for k in 1:div(D,2)
        cache(int,ST).q̃[k] = (sol.q[k] + cache(int,ST).q[k]) / 2
        cache(int,ST).q̃[div(D,2)+k] = x[D+div(D,2)+k]

        cache(int,ST).v[k]          = x[D+k]
        cache(int,ST).v[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    equations(int).f(cache(int,ST).f̃, t̃, cache(int,ST).q̃, cache(int,ST).v, params)
 
    # compute Θ = ϑ(q,v)
    equations(int).ϑ(cache(int,ST).p̃, t̃, cache(int,ST).q̃, cache(int,ST).v, params)
    equations(int).ϑ(cache(int,ST).p, t, cache(int,ST).q, cache(int,ST).v, params)
end


function residual!(b::Vector{ST}, sol, params, int::GeometricIntegrator{<:CMDVI}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # compute b
    b[1:D] .= cache(int,ST).p̃ .- sol.p .- timestep(int) .* cache(int,ST).f̃ ./ 2
    
    for k in 1:div(D,2)
        b[D+k]          = cache(int,ST).q[k] - sol.q[k] - timestep(int) * cache(int,ST).v[k]
        b[D+div(D,2)+k] = cache(int,ST).p[k] - sol.p[k] - timestep(int) * cache(int,ST).f̃[k]
    end
end


function update!(sol, params, int::GeometricIntegrator{<:CMDVI}, DT)
    # compute final update
    sol.q .= cache(int, DT).q
    sol.p .= cache(int, DT).p
end

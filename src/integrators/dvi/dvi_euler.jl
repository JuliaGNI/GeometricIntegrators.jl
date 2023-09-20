
abstract type DVIEuler <: DVIMethod end

"Symplectic Euler-A Degenerate Variational Integrator."
struct DVIA <: DVIEuler end

"Symplectic Euler-B Degenerate Variational Integrator."
struct DVIB <: DVIEuler end

order(::Union{DVIA, Type{DVIA}}) = 1
order(::Union{DVIB, Type{DVIB}}) = 1

issymmetric(::Union{DVIA, Type{<:DVIA}}) = false
issymmetric(::Union{DVIB, Type{<:DVIB}}) = false


function Base.show(io::IO, int::GeometricIntegrator{<:DVIA})
    print(io, "\nDegenerate Variational Integrator (Euler-A) with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end

function Base.show(io::IO, int::GeometricIntegrator{<:DVIB})
    print(io, "\nDegenerate Variational Integrator (Euler-B) with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
end


function initial_guess!(int::GeometricIntegrator{<:DVIEuler})
    # set some local variables for convenience
    local D = ndims(int)
    local x = nlsolution(int)

    # compute initial guess for solution
    initialguess!(solstep(int).t, cache(int).q, cache(int).p, cache(int).v, cache(int).f, solstep(int), problem(int), iguess(int))

    # copy q to nonlinear solution vector
    x[1:D] .= cache(int).q

    # copy v to nonlinear solution vector
    for k in 1:div(D,2)
        x[D+k] = cache(int).v[k]
        x[D+div(D,2)+k] = solstep(int).v[k]
    end
end


function components!(x::Vector{ST}, int::GeometricIntegrator{<:DVIEuler}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # copy x to q
    cache(int,ST).q .= x[1:D]

    # copy x to v and v̄
    for k in 1:div(D,2)
        cache(int,ST).v[k] = x[D+k]
        cache(int,ST).v̄[k] = x[D+div(D,2)+k]
        cache(int,ST).v[div(D,2)+k] = 0
        cache(int,ST).v̄[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    equations(int).f(cache(int,ST).f, solstep(int).t, cache(int,ST).q, cache(int,ST).v, parameters(solstep(int)))
    equations(int).f(cache(int,ST).f̄, solstep(int).t̄, cache(int).q̄, cache(int,ST).v̄, parameters(solstep(int)))
 
    # compute Θ = ϑ(q,v)
    equations(int).ϑ(cache(int,ST).p, solstep(int).t, cache(int,ST).q, cache(int,ST).v, parameters(solstep(int)))
    # equations(int).ϑ(cache(int,ST).θ̄, solstep(int).t̄, solstep(int).q̄, cache(int,ST).v̄, parameters(solstep(int)))
end


function residual!(b::Vector{ST}, x::Vector{ST}, int::GeometricIntegrator{<:DVIA}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute b
    for k in 1:div(D,2)
        b[k]   = cache(int,ST).p[k] - cache(int).p̄[k] - timestep(int) * cache(int,ST).f̄[k]
        b[D+k] = cache(int,ST).q[k] - cache(int).q̄[k] - timestep(int) * cache(int,ST).v[k]
    end

    for k in div(D,2)+1:D
        b[k]   = timestep(int) * cache(int,ST).f̄[k]
        b[D+k] = timestep(int) * cache(int,ST).f[k]
    end
end


function residual!(b::Vector{ST}, x::Vector{ST}, int::GeometricIntegrator{<:DVIB}) where {ST}
    # set some local variables for convenience
    local D = ndims(int)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute b
    for k in 1:div(D,2)
        b[k]   = cache(int,ST).p[k] - cache(int).p̄[k] - timestep(int) * cache(int,ST).f[k]
        b[D+k] = cache(int,ST).q[k] - cache(int).q̄[k] - timestep(int) * cache(int,ST).v̄[k]
    end

    for k in div(D,2)+1:D
        b[k]   = timestep(int) * cache(int,ST).f[k]
        b[D+k] = timestep(int) * cache(int,ST).f̄[k]
    end
end


function update!(x::AbstractVector{DT}, int::GeometricIntegrator{<:DVIEuler}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    solstep(int).q .= cache(int, DT).q
    solstep(int).v .= cache(int, DT).v
    solstep(int).p .= cache(int, DT).p
    solstep(int).f .= cache(int, DT).f
end

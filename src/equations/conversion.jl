
@define _create_pode_argument_views begin
    n = div(length(eachindex(x)), 2)
    q = @view x[eachindex(x)[  1:n ]]
    p = @view x[eachindex(x)[n+1:2n]]
    q̇ = @view ẋ[eachindex(ẋ)[  1:n ]]
    ṗ = @view ẋ[eachindex(ẋ)[n+1:2n]]
end


convert_periodicity(::Union{Type{ODE}, Type{SODE}}, equ::Union{PODE, HODE}) = vcat(periodicity(equ), zero(periodicity(equ)))


function Base.convert(::Type{ODE}, equ::Union{PODE{DT,TT,AT}, HODE{DT,TT,AT}}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = [vcat(x...) for x in zip(equ.q₀, equ.p₀)]

    # extend periodicity
    ode_periodicity = convert_periodicity(ODE, equ)

    if hasparameters(equ)
        v = (t, x, ẋ, params) -> begin
            @_create_pode_argument_views
            equ.v(t, q, p, q̇, params)
            equ.f(t, q, p, ṗ, params)
        end
    else
        v = (t, x, ẋ) -> begin
            @_create_pode_argument_views
            equ.v(t, q, p, q̇)
            equ.f(t, q, p, ṗ)
        end
    end

    ODE(v, equ.t₀, x₀; parameters=equ.parameters, periodicity=ode_periodicity)
    # TODO: Convert invariants and pass to ODE
    # TODO: For HODE append (h=equ.h,) to invariants
end

function Base.convert(::Type{SODE}, equ::Union{PODE{DT,TT,AT}, HODE{DT,TT,AT}}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = [vcat(x...) for x in zip(equ.q₀, equ.p₀)]

    # extend periodicity
    ode_periodicity = convert_periodicity(SODE, equ)

    if hasparameters(equ)
        v₁ = (t, x, ẋ, params) -> begin
            @_create_pode_argument_views
            equ.v(t, q, p, q̇, params)
        end
        v₂ = (t, x, ẋ, params) -> begin
            @_create_pode_argument_views
            equ.f(t, q, p, ṗ, params)
        end
    else
        v₁ = (t, x, ẋ) -> begin
            @_create_pode_argument_views
            equ.v(t, q, p, q̇)
        end
        v₂ = (t, x, ẋ) -> begin
            @_create_pode_argument_views
            equ.f(t, q, p, ṗ)
        end
    end

    SODE((v₁, v₂), equ.t₀, x₀; parameters=equ.parameters, periodicity=ode_periodicity)
    # TODO: Convert invariants and pass to SODE
end

function Base.convert(::Type{PODE}, equ::HODE)
    PODE(equ.v, equ.f, equ.t₀, equ.q₀, equ.p₀;
         invariants=equ.invariants, parameters=equ.parameters, periodicity=equ.periodicity)
    # TODO: Append (h=equ.h,) to invariants
end

function Base.convert(::Type{IODE}, equ::LODE)
    IODE(equ.ϑ, equ.f, equ.g, equ.t₀, equ.q₀, equ.p₀, equ.λ₀;
         v̄=equ.v̄, f̄=equ.f̄, invariants=equ.invariants, parameters=equ.parameters, periodicity=equ.periodicity)
end


function get_invariants(equ::Union{ODE,SODE,DAE})
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t,q) -> inv(t, q, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end

function get_invariants(equ::Union{IODE,LODE,IDAE,LDAE})
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t,q,v) -> inv(t, q, v, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end

function get_invariants(equ::Union{PODE,HODE,PDAE,PDAE,SPDAE})
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t,q,p) -> inv(t, q, p, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end

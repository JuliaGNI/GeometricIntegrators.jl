
@define _create_pode_argument_views begin
    n = div(length(eachindex(x)), 2)
    q = @view x[eachindex(x)[  1:n ]]
    p = @view x[eachindex(x)[n+1:2n]]
    q̇ = @view ẋ[eachindex(ẋ)[  1:n ]]
    ṗ = @view ẋ[eachindex(ẋ)[n+1:2n]]
end


function Base.convert(::Type{ODE}, equ::Union{PODE{DT,TT,AT}, HODE{DT,TT,AT}}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = [vcat(x...) for x in zip(equ.q₀, equ.p₀)]

    # extend periodicity
    periodicity = vcat(equ.periodicity, zero(equ.periodicity))

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

    ODE(v, equ.t₀, x₀; h=equ.h, parameters=equ.parameters, periodicity=periodicity)
end

function Base.convert(::Type{SODE}, equ::Union{PODE{DT,TT,AT}, HODE{DT,TT,AT}}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = [vcat(x...) for x in zip(equ.q₀, equ.p₀)]

    # extend periodicity
    periodicity = vcat(equ.periodicity, zero(equ.periodicity))

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

    SODE((v₁, v₂), equ.t₀, x₀; parameters=equ.parameters, periodicity=periodicity)
end

function Base.convert(::Type{PODE}, equ::HODE)
    PODE(equ.v, equ.f, equ.t₀, equ.q₀, equ.p₀; h=equ.h, parameters=equ.parameters, periodicity=equ.periodicity)
end

function Base.convert(::Type{IODE}, equ::LODE)
    IODE(equ.ϑ, equ.f, equ.g, equ.t₀, equ.q₀, equ.p₀, equ.λ₀; v̄=equ.v̄, f̄=equ.f̄, h=equ.h, parameters=equ.parameters, periodicity=equ.periodicity)
end

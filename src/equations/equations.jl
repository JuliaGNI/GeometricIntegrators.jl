
abstract type Equation{dType <: Number, tType <: Real} end

abstract type AbstractEquationODE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPODE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationDAE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPDAE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationSDE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPSDE{dType, tType} <: Equation{dType, tType} end


Base.eltype(::Equation{DT,TT}) where {DT,TT} = DT
timetype(::Equation{DT,TT}) where {DT,TT} = TT


function_v_dummy(t, q, p, v) = nothing

function get_λ₀(q₀::AbstractVector{DT}, λ₀::AbstractVector{DT}) where {DT <: Number}
    zero(λ₀)
end

function get_λ₀(q₀::AbstractVector{DT}, λ₀::AbstractVector{AT}) where {DT, AT <: AbstractArray{DT}}
    zero(λ₀[begin])
end

function get_λ₀(q₀::AbstractVector{AT}, λ₀::AbstractVector{DT}) where {DT, AT <: AbstractArray{DT}}
    [zero(λ₀) for i in eachindex(q₀)]
end

function get_λ₀(q₀::AbstractVector{AT}, λ₀::AbstractVector{AT}) where {DT, AT <: AbstractArray{DT}}
    [zero(λ₀[begin]) for i in eachindex(q₀)]
end

Base.ndims(equ::Equation) = error("ndims() not implemented for ", typeof(equ), ".")

periodicity(equ::Equation) = error("periodicity() not implemented for ", typeof(equ), ".")

get_function_tuple(equ::Equation) = error("get_function_tuple() not implemented for ", typeof(equ), ".")

# TODO Add functions and vectors for invariants (momentum maps, energy, ...).

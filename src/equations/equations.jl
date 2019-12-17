
abstract type Equation{dType <: Number, tType <: Real} end

abstract type AbstractEquationODE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPODE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationDAE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPDAE{dType, tType} <: Equation{dType, tType} end


function function_v_dummy(t, q, p, v)
    nothing
end

function get_λ₀(q₀::AbstractArray{DT,1}, λ₀) where {DT}
    zeros(DT, size(λ₀,1))
end

function get_λ₀(q₀::AbstractArray{DT,2}, λ₀) where {DT}
    zeros(DT, size(λ₀,1), size(q₀,2))
end

Base.ndims(equ::Equation) = error("ndims() not implemented for ", typeof(equ), ".")

periodicity(equ::Equation) = error("periodicity() not implemented for ", typeof(equ), ".")


# TODO Add functions and vectors for invariants (momentum maps, energy, ...).

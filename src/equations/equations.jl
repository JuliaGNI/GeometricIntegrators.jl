
abstract type Equation{dType <: Number, tType <: Real} end

function function_v_dummy(t, q, p, v)
    nothing
end

Base.ndims(equ::Equation) = error("ndims() not implemented for ", typeof(equ))

# TODO Add functions and vectors for invariants (momentum maps, energy, ...).

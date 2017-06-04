
abstract type Equation{dType <: Number, tType <: Real} end

function function_v_dummy(t, q, p, v)
    nothing
end

# TODO Add functions and vectors for invariants (momentum maps, energy, ...).

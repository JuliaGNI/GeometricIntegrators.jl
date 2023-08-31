"""
Explicit Euler Method.

$(reference(Val(:ExplicitEuler)))
"""
struct ExplicitEuler <: ODEMethod end

isexplicit(method::ExplicitEuler) = true
isimplicit(method::ExplicitEuler) = false
issymmetric(method::ExplicitEuler) = false
issymplectic(method::ExplicitEuler) = false

const ExplicitEulerIntegrator{DT,TT} = GeometricIntegrator{<:ExplicitEuler, <:AbstractProblemODE{DT,TT}}


function integrate_step!(int::ExplicitEulerIntegrator)
    # compute vector field
    equations(int)[:v](solstep(int).v, solstep(int).t, solstep(int).q)

    # compute update
    update!(solstep(int), solstep(int).v, timestep(int))
end

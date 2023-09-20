"""
Explicit Euler Method.

$(reference(Val(:ExplicitEuler)))
"""
struct ExplicitEuler <: ODEMethod end

isexplicit(method::ExplicitEuler) = true
isimplicit(method::ExplicitEuler) = false
issymmetric(method::ExplicitEuler) = false
issymplectic(method::ExplicitEuler) = false


function integrate_step!(int::GeometricIntegrator{<:ExplicitEuler, <:AbstractProblemODE})
    # compute vector field
    equations(int).v(solstep(int).v, solstep(int).t, solstep(int).q, parameters(solstep(int)))

    # compute update
    update!(solstep(int), solstep(int).v, timestep(int))
end

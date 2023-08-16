"""
Explicit Euler Method.

$(reference(Val(:ExplicitEuler)))
"""
struct ExplicitEuler <: ODEMethod end

Methods.isexplicit(method::ExplicitEuler) = true
Methods.isimplicit(method::ExplicitEuler) = false
Methods.issymmetric(method::ExplicitEuler) = false
Methods.issymplectic(method::ExplicitEuler) = false


const IntegratorExplicitEuler{DT,TT} = GeometricIntegrator{<:Union{ODEProblem{DT,TT}, DAEProblem{DT,TT}, SubstepProblem{DT,TT}}, <:ExplicitEuler}

function integrate_step!(int::IntegratorExplicitEuler)

    # compute vector field
    equations(int)[:v](solstep(int).v, solstep(int).t, solstep(int).q)

    # compute update
    update!(solstep(int), solstep(int).v, timestep(int))
end

"""
Explicit Euler Method.

$(reference(Val(:ExplicitEuler)))
"""
struct ExplicitEuler <: ODEMethod end

const IntegratorExplicitEuler{DT,TT} = Integrator{<:Union{ODEProblem{DT,TT}, DAEProblem{DT,TT}, SubstepProblem{DT,TT}}, <:ExplicitEuler}

function integrate_step!(int::IntegratorExplicitEuler)

    # compute vector field
    equations(int)[:v](solstep(int).v, solstep(int).t̄, solstep(int).q)

    # compute update
    update!(solstep(int), solstep(int).v, timestep(int))
end

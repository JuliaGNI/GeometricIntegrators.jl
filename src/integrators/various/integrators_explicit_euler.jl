
"Explicit Euler integrator."
struct IntegratorExplicitEuler{DT, TT, AT, D, ET <: NamedTuple} <: DeterministicIntegrator{DT,TT}
    equs::ET
    Δt::TT
    v::AT

    function IntegratorExplicitEuler{DT,AT,D}(equations::ET, Δt::TT) where {DT, TT, AT, D, ET}
        # create cache for vector field
        v = AT(zeros(DT, D))

        # create integrator
        new{DT,TT,AT,D,ET}(equations, Δt, v)
    end

    function IntegratorExplicitEuler{DT,AT,D}(v::Function, Δt::TT; kwargs...) where {DT,TT,AT,D}
        IntegratorExplicitEuler{DT,AT,D}(NamedTuple{(:v,)}((v,)), Δt; kwargs...)
    end

    function IntegratorExplicitEuler(equation::ODE{DT,TT,AT}, Δt::TT; kwargs...) where {DT,TT,AT}
        IntegratorExplicitEuler{DT, AT, axes(equation)}(get_function_tuple(equation), Δt; kwargs...)
    end
end

equation(int::IntegratorExplicitEuler, i::Symbol) = int.equs[i]
equations(int::IntegratorExplicitEuler) = int.equs
timestep(int::IntegratorExplicitEuler) = int.Δt


function integrate_step!(int::IntegratorExplicitEuler{DT,TT,AT}, sol::AtomicSolutionODE{DT,TT,AT}) where {DT,TT,AT}
    # reset atomic solution
    reset!(sol, timestep(int))

    # compute vector field
    equations(int)[:v](sol.t̄, sol.q, int.v)

    # compute update
    sol.q .+= timestep(int) .* int.v
end

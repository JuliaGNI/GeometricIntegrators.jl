
"Exact solution of an ODE."
struct IntegratorExactODE{DT, TT, D, QT <: Function} <: ODEIntegrator{DT,TT}
    q::QT
    Δt::TT

    function IntegratorExactODE{DT,D}(solution::solType, Δt::TT) where {DT, TT, D, solType <: Function}
        new{DT,TT,D,solType}(solution, Δt)
    end
end


@inline Base.ndims(::IntegratorExactODE{DT,TT,D}) where {DT,TT,D} = D
@inline GeometricBase.timestep(int::IntegratorExactODE) = int.Δt


function integrate_step!(int::IntegratorExactODE{DT,TT}, sol::SolutionStepODE{DT,TT}) where {DT,TT}
    # reset atomic solution
    reset!(sol)

    # compute new solution
    int.q(sol.q, sol.t̄ + timestep(int), sol.q̄, sol.t̄)
end

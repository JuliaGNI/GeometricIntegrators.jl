
"Exact solution of an ODE."
struct IntegratorExactODE{DT, TT, D, QT <: Function} <: ODEIntegrator{DT,TT}
    q::QT
    Δt::TT

    function IntegratorExactODE{DT,D}(solution::solType, Δt::TT) where {DT, TT, D, solType <: Function}
        new{DT,TT,D,solType}(solution, Δt)
    end
end


@inline Base.ndims(::IntegratorExactODE{DT,TT,D}) where {DT,TT,D} = D
timestep(int::IntegratorExactODE) = int.Δt


function integrate_step!(int::IntegratorExactODE{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    # reset atomic solution
    reset!(sol, timestep(int))

    # compute new solution
    int.q(sol.t̅ + timestep(int), sol.q̅, sol.q, timestep(int))
end

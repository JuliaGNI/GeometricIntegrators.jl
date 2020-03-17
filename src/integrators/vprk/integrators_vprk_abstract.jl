
abstract type AbstractIntegratorVPRK{DT,TT,D,S} <: DeterministicIntegrator{DT,TT} end
abstract type AbstractIntegratorVPRKwProjection{DT,TT,D,S} <: AbstractIntegratorVPRK{DT,TT,D,S} end

@inline parameters(integrator::AbstractIntegratorVPRK) = integrator.params

@inline nstages(integrator::AbstractIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = S
@inline Base.ndims(integrator::AbstractIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = D

@inline eachstage(integrator::AbstractIntegratorVPRK) = 1:nstages(integrator)
@inline eachdim(integrator::AbstractIntegratorVPRK) = 1:ndims(integrator)


function update_solution!(int::AbstractIntegratorVPRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    update_solution!(sol.q, sol.q̃, int.cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))
end

function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT}, R::Vector{TT}) where {DT,TT}
    update_solution!(sol.q, sol.q̃, int.cache.U, R, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.G, R, timestep(int))
end

function project_solution!(int::AbstractIntegratorVPRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT}, RU::Vector{TT}, RG::Vector{TT}) where {DT,TT}
    update_solution!(sol.q, sol.q̃, int.cache.U, RU, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.G, RG, timestep(int))
end

function Integrators.initialize!(int::AbstractIntegratorVPRK{DT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    sol.t̅ = sol.t - timestep(int)

    equation(int, :v)(sol.t, sol.q, sol.p, sol.v)
    equation(int, :f)(sol.t, sol.q, sol.p, sol.f)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̅, sol.q̅, sol.p̅, sol.v̅, sol.f̅)
end

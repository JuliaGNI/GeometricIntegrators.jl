
abstract type IntegratorCache{DT,D} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


create_integrator_cache(integrator::Integrator) = error("create_integrator_cache()! not implemented for ", typeof(integrator))

initialize!(::Integrator, ::IntegratorCache) = nothing

integrate_step!(integrator::Integrator, ::IntegratorCache) = error("integrate_step()! not implemented for ", typeof(integrator))

function copy_solution!(sol::Solution, cache::IntegratorCache, n, m)
    copy_solution!(sol, get_solution(cache)..., n, m)
end

function CommonFunctions.get_solution(cache::ODEIntegratorCache)
    (cache.t, cache.q)
end

function CommonFunctions.get_solution(cache::Union{PODEIntegratorCache, IODEIntegratorCache})
    (cache.t, cache.q, cache.p)
end

function CommonFunctions.get_solution(cache::Union{PDAEIntegratorCache, IDAEIntegratorCache})
    (cache.t, cache.q, cache.p, cache.λ)
end

function cut_periodic_solution!(cache::IntegratorCache, periodicity::Vector)
    cut_periodic_solution!(cache.q, periodicity, cache.s̃)
    cache.q .+= cache.s̃
    cache.q̅ .+= cache.s̃
end

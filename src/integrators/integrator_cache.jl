
abstract type IntegratorCache{DT,D} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


create_integrator_cache(integrator::Integrator) = error("create_integrator_cache()! not implemented for ", typeof(integrator))

initialize!(::Integrator, ::IntegratorCache) = nothing


function copy_solution!(sol::Solution, cache::IntegratorCache, n, m)
    copy_solution!(sol, get_solution(cache)..., n, m)
end

function CommonFunctions.get_solution(cache::ODEIntegratorCache)
    (cache.t, cache.q)
end

function CommonFunctions.get_solution(cache::PODEIntegratorCache)
    (cache.t, cache.q, cache.p)
end

function cut_periodic_solution!(cache::Union{ODEIntegratorCache, PODEIntegratorCache}, periodicity::Vector)
    cut_periodic_solution!(cache.q, periodicity, cache.s̃)
    cache.q .+= cache.s̃
    cache.q̅ .+= cache.s̃
end

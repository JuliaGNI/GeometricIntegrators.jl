
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


function CommonFunctions.reset!(cache::Union{ODEIntegratorCache, DAEIntegratorCache}, Δt)
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.v̅ .= cache.v
    cache.t += Δt
    cache.n += 1
end

function CommonFunctions.reset!(cache::Union{PODEIntegratorCache, IODEIntegratorCache}, Δt)
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.p̅ .= cache.p
    cache.v̅ .= cache.v
    cache.f̅ .= cache.f
    cache.t += Δt
    cache.n += 1
end

function CommonFunctions.reset!(cache::Union{PDAEIntegratorCache, IDAEIntegratorCache}, Δt)
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.p̅ .= cache.p
    cache.λ̅ .= cache.λ
    cache.v̅ .= cache.v
    cache.f̅ .= cache.f
    cache.u̅ .= cache.u
    cache.g̅ .= cache.g
    cache.t += Δt
    cache.n += 1
end


function CommonFunctions.set_solution!(cache::ODEIntegratorCache, sol, n=0)
    t, q = sol
    cache.n  = n
    cache.t  = t
    cache.q .= q
    cache.v .= 0
end


function CommonFunctions.set_solution!(cache::Union{PODEIntegratorCache, IODEIntegratorCache}, sol, n=0)
    t, q, p = sol
    cache.n  = n
    cache.t  = t
    cache.q .= q
    cache.p .= p
    cache.v .= 0
    cache.f .= 0
end

function CommonFunctions.set_solution!(cache::Union{PDAEIntegratorCache, IDAEIntegratorCache}, sol, n=0)
    t, q, p, λ = sol
    cache.n  = n
    cache.t  = t
    cache.q .= q
    cache.p .= p
    cache.λ .= λ
    cache.v .= 0
    cache.f .= 0
end


function CommonFunctions.get_solution(cache::ODEIntegratorCache)
    (cache.t, cache.q)
end

function CommonFunctions.get_solution(cache::DAEIntegratorCache)
    (cache.t, cache.q, cache.λ)
end

function CommonFunctions.get_solution(cache::Union{PODEIntegratorCache, IODEIntegratorCache})
    (cache.t, cache.q, cache.p)
end

function CommonFunctions.get_solution(cache::Union{PDAEIntegratorCache, IDAEIntegratorCache})
    (cache.t, cache.q, cache.p, cache.λ)
end


function copy_solution!(sol::Solution, cache::IntegratorCache, n, m)
    copy_solution!(sol, get_solution(cache)..., n, m)
end


function cut_periodic_solution!(cache::IntegratorCache, periodicity::Vector)
    cut_periodic_solution!(cache.q, periodicity, cache.s̃)
    cache.q .+= cache.s̃
    cache.q̅ .+= cache.s̃
end

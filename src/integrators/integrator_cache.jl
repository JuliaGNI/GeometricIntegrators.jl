
abstract type IntegratorCache{DT,D} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


IntegratorCache(params::Parameters{DT}) where {DT} = IntegratorCache{DT}(params)
IntegratorCache{ST}(params::Parameters) where {ST} = error("IntegratorCache{ST}(params) not implemented for ", typeof(params))

IntegratorCache(integrator::Integrator{DT}) where {DT} = IntegratorCache{DT}(integrator)
IntegratorCache{ST}(integrator::Integrator) where {ST} = error("IntegratorCache{ST}(int)! not implemented for ", typeof(integrator))

CacheType(T, params::Parameters) = error("CacheType(T, params) not implemented for ", typeof(params))


struct CacheDict{PT}
    params::PT
    caches::Dict{UInt64, IntegratorCache}

    function CacheDict(params::Parameters)
        new{typeof(params)}(params, Dict{UInt64, IntegratorCache}())
    end
end

@inline function Base.getindex(c::CacheDict, ST::DataType)
    key = hash(Threads.threadid(), hash(ST))
    if haskey(c.caches, key)
        c.caches[key]
    else
        c.caches[key] = IntegratorCache{ST}(c.params)
    end::CacheType(ST, c.params)
end


abstract type Cache{DT} end

abstract type IntegratorCache{DT,D} <: Cache{DT} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


IntegratorCache(params::Parameters{DT}) where {DT} = IntegratorCache{DT}(params)
IntegratorCache{ST}(params::Parameters) where {ST} = error("IntegratorCache{ST}(params) not implemented for ", typeof(params))

IntegratorCache(integrator::AbstractIntegrator{DT}) where {DT} = IntegratorCache{DT}(integrator)
IntegratorCache{ST}(integrator::AbstractIntegrator) where {ST} = error("IntegratorCache{ST}(int) not implemented for ", typeof(integrator))

CacheType(T, params::Parameters) = error("CacheType(T, params) not implemented for ", typeof(params))

struct OldCacheDict{PT}
    params::PT
    caches::Dict{UInt64, IntegratorCache}

    function OldCacheDict(params::Parameters)
        new{typeof(params)}(params, Dict{UInt64, IntegratorCache}())
    end
end

@inline function Base.getindex(c::OldCacheDict, ST::DataType)
    key = hash(Threads.threadid(), hash(ST))
    if haskey(c.caches, key)
        c.caches[key]
    else
        c.caches[key] = IntegratorCache{ST}(c.params)
    end::CacheType(ST, c.params)
end


IntegratorCache(problem::GeometricProblem, method::GeometricMethod) = IntegratorCache{datatype(problem)}(problem, method)
IntegratorCache{ST}(::GeometricProblem, ::GeometricMethod) where {ST} = nothing

CacheType(T, problem::GeometricProblem, method::GeometricMethod) = error("CacheType(T, params) not implemented for ", typeof(problem), " and ", typeof(method))


struct CacheDict{PT,MT}
    problem::PT
    method::MT
    caches::Dict{UInt64, IntegratorCache}

    function CacheDict(prob::GeometricProblem, method::GeometricMethod)
        new{typeof(prob), typeof(method)}(prob, method, Dict{UInt64, IntegratorCache}())
    end
end

@inline function Base.getindex(c::CacheDict, ST::DataType)
    key = hash(Threads.threadid(), hash(ST))
    if haskey(c.caches, key)
        c.caches[key]
    else
        c.caches[key] = Cache{ST}(c.problem, c.method)
    end::CacheType(ST, c.problem, c.method)
end

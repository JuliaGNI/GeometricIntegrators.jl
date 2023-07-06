
abstract type Cache{DT} end

abstract type IntegratorCache{DT,D} <: Cache{DT} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


reset!(::ODEIntegratorCache, t, q, λ = missing) = nothing


IntegratorCache(problem::AbstractProblem, method::GeometricMethod) = IntegratorCache{datatype(problem)}(problem, method)
IntegratorCache{ST}(::AbstractProblem, ::GeometricMethod) where {ST} = nothing

CacheType(T, problem::AbstractProblem, method::GeometricMethod) = error("CacheType(T, params) not implemented for ", typeof(problem), " and ", typeof(method))


struct CacheDict{PT,MT}
    problem::PT
    method::MT
    caches::Dict{UInt64, IntegratorCache}

    function CacheDict(prob::AbstractProblem, method::GeometricMethod)
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

function updateall!(c::CacheDict, data)
    for (key, cache) in c.caches
        update!(c.caches[key], data...)
    end
end

cache(c::CacheDict) = c[datatype(c.problem)]
nlsolution(c::CacheDict) = nlsolution(cache(c))

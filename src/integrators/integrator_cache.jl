
abstract type IntegratorCache{DT,D,S} end

abstract type ODEIntegratorCache{DT,D,S} <: IntegratorCache{DT,D,S} end
abstract type DAEIntegratorCache{DT,D,S} <: IntegratorCache{DT,D,S} end
abstract type IODEIntegratorCache{DT,D,S} <: IntegratorCache{DT,D,S} end
abstract type IDAEIntegratorCache{DT,D,S} <: IntegratorCache{DT,D,S} end
abstract type PODEIntegratorCache{DT,D,S} <: IntegratorCache{DT,D,S} end
abstract type PDAEIntegratorCache{DT,D,S} <: IntegratorCache{DT,D,S} end


create_integrator_cache(integrator::Integrator) = error("create_integrator_cache()! not implemented for ", typeof(integrator))


function copy_solution!(sol::Solution, cache::IntegratorCache, n, m)
    copy_solution!(sol, get_solution(cache)..., n, m)
end


# function initialize!(int::IntegratorERK, sol::SolutionODE, m::Int)
#     @assert m ≥ 1
#     @assert m ≤ sol.ni
#
#     # copy initial conditions from solution
#     get_initial_conditions!(sol, int.q, m)
# end

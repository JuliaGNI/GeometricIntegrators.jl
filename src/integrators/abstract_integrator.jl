
abstract type Integrator{dType, tType} end

# function CommonFunctions.name(int::Integrator)
#     warn(string(typeof(int)) * ".name() Not implemented!")
# end

equation(integrator::Integrator) = integrator.equation
timestep(integrator::Integrator) = integrator.Δt
integrate(integrator::Integrator) = error("integrate() not implemented for ", typeof(integrator))
integrate!(integrator::Integrator) = error("integrate()! not implemented for ", typeof(integrator))


abstract type NonlinearFunctionParameters{DT,TT} end

function_stages!(x::Vector{DT}, b::Vector{DT}, params::PT) where {DT, TT, PT <: NonlinearFunctionParameters{DT,TT}} = error("function_stages!() not implemented for ", PT)

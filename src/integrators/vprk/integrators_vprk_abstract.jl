
# abstract type AbstractIntegratorVPRK{DT,TT,D,S} <: AbstractIntegratorIRK{DT,TT} end
# abstract type AbstractIntegratorVPRKwProjection{DT,TT,D,S} <: AbstractIntegratorVPRK{DT,TT,D,S} end

# @inline parameters(integrator::AbstractIntegratorVPRK) = integrator.params

# @inline nstages(::AbstractIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = S
# @inline Base.ndims(::AbstractIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = D

# @inline eachstage(integrator::AbstractIntegratorVPRK) = 1:nstages(integrator)
# @inline eachdim(integrator::AbstractIntegratorVPRK) = 1:ndims(integrator)

# function Base.show(io::IO, int::AbstractIntegratorVPRK)
#     print(io, "\n$(description(int)) with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(method(int)))\n")
#     print(io, "   $(string(method(int).q))")
#     print(io, "   $(string(method(int).p))")
#     # print(io, reference(method(int)))
# end

# function Base.show(io::IO, int::AbstractIntegratorVPRKwProjection)
#     print(io, "\n$(description(int)) and:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(method(int)))\n")
#     print(io, "   $(string(method(int).q))")
#     print(io, "   $(string(method(int).p))")
#     # print(io, reference(method(int)))
# end




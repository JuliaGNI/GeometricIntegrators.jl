
# abstract type GeometricIntegratorVPRK{DT,TT,D,S} <: AbstractIntegratorIRK{DT,TT} end
# abstract type GeometricIntegratorVPRKwProjection{DT,TT,D,S} <: GeometricIntegratorVPRK{DT,TT,D,S} end

# @inline parameters(integrator::GeometricIntegratorVPRK) = integrator.params

# @inline nstages(::GeometricIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = S
# @inline Base.ndims(::GeometricIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = D

# @inline eachstage(integrator::GeometricIntegratorVPRK) = 1:nstages(integrator)
# @inline eachdim(integrator::GeometricIntegratorVPRK) = 1:ndims(integrator)

# function Base.show(io::IO, int::GeometricIntegratorVPRK)
#     print(io, "\n$(description(int)) with:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(method(int)))\n")
#     print(io, "   $(string(method(int).q))")
#     print(io, "   $(string(method(int).p))")
#     # print(io, reference(method(int)))
# end

# function Base.show(io::IO, int::GeometricIntegratorVPRKwProjection)
#     print(io, "\n$(description(int)) and:\n")
#     print(io, "   Timestep: $(timestep(int))\n")
#     print(io, "   Tableau:  $(description(method(int)))\n")
#     print(io, "   $(string(method(int).q))")
#     print(io, "   $(string(method(int).p))")
#     # print(io, reference(method(int)))
# end




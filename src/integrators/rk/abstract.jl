
const StageVector{T} = Vector{<:AbstractVector{T}}

GeometricBase.tableau(int::GeometricIntegrator{<:EquationProblem, <:RKMethod}) = tableau(method(int))

initmethod(method::RKMethod) = RK(method)
initmethod(method::PRKMethod) = PRK(method)

GeometricIntegrator(problem::Union{PODEProblem,HODEProblem}, method::RKMethod, args...; kwargs...) = GeometricIntegrator(problem, PRK(method), args...; kwargs...)

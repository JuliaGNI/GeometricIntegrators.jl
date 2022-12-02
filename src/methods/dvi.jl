
abstract type DVIMethod <: LODEMethod end


"Symplectic Euler-A Degenerate Variational Integrator."
struct DVIA <: DVIMethod end

"Symplectic Euler-B Degenerate Variational Integrator."
struct DVIB <: DVIMethod end

"Midpoint Degenerate Variational Integrator."
struct CMDVI <: DVIMethod end

"Trapezoidal Degenerate Variational Integrator."
struct CTDVI <: DVIMethod end

order(::Union{DVIA, Type{DVIA}}) = 1
order(::Union{DVIB, Type{DVIB}}) = 1
order(::Union{CMDVI, Type{CMDVI}}) = 2
order(::Union{CTDVI, Type{CTDVI}}) = 2

issymmetric(::Union{DVIA, Type{<:DVIA}}) = false
issymmetric(::Union{DVIB, Type{<:DVIB}}) = false
issymmetric(::Union{CMDVI, Type{<:CMDVI}}) = true
issymmetric(::Union{CTDVI, Type{<:CTDVI}}) = true

issymplectic(::Union{DVIMethod, Type{<:DVIMethod}}) = true

Integrators.Integrator(problem::LODEProblem, method::DVIA; kwargs...) = IntegratorDVIA(problem; kwargs...)
Integrators.Integrator(problem::LODEProblem, method::DVIB; kwargs...) = IntegratorDVIB(problem; kwargs...)
Integrators.Integrator(problem::LODEProblem, method::CMDVI; kwargs...) = IntegratorCMDVI(problem; kwargs...)
Integrators.Integrator(problem::LODEProblem, method::CTDVI; kwargs...) = IntegratorCTDVI(problem; kwargs...)

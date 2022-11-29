
abstract type DVIMethod <: LODEMethod end


"Symplectic Euler-A Degenerate Variational Integrator."
struct DVIA <: DVIMethod end

"Symplectic Euler-B Degenerate Variational Integrator."
struct DVIB <: DVIMethod end

"Midpoint Degenerate Variational Integrator."
struct CMDVI <: DVIMethod end

"Trapezoidal Degenerate Variational Integrator."
struct CTDVI <: DVIMethod end


Integrators.Integrator(problem::LODEProblem, method::DVIA; kwargs...) = IntegratorDVIA(problem; kwargs...)
Integrators.Integrator(problem::LODEProblem, method::DVIB; kwargs...) = IntegratorDVIB(problem; kwargs...)
Integrators.Integrator(problem::LODEProblem, method::CMDVI; kwargs...) = IntegratorCMDVI(problem; kwargs...)
Integrators.Integrator(problem::LODEProblem, method::CTDVI; kwargs...) = IntegratorCTDVI(problem; kwargs...)

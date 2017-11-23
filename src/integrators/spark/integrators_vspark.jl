"Holds the tableau of an Variational Specialised Partitioned Additive Runge-Kutta method."
struct TableauVSPARK{T} <: AbstractTableau{T}
end

"Parameters for right-hand side function of Variational Specialised Partitioned Additive Runge-Kutta methods."
mutable struct NonlinearFunctionParametersVSPARK{DT,TT,VT,FT,ϕT,ψT} <: NonlinearFunctionParameters{DT,TT}
end


"Variational Specialised Partitioned Additive Runge-Kutta integrator."
immutable IntegratorVSPARK{DT, TT, VT, FT, ϕT, ψT, SPT, ST, IT} <: Integrator{DT, TT}
end

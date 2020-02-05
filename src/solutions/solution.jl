
using HDF5

abstract type Solution{dType, tType, N} end
abstract type DeterministicSolution{dType, tType, N} <: Solution{dType, tType, N} end
abstract type StochasticSolution{dType, tType, NQ, NW} <: Solution{dType, tType, NQ} end

timesteps(sol::Solution) = error("time() not implemented for ", typeof(sol))
hdf5(sol::Solution)   = error("hdf5() not implemented for ", typeof(sol))
ntime(sol::Solution)  = error("ntime() not implemented for ", typeof(sol))
nsave(sol::Solution)  = error("nsave() not implemented for ", typeof(sol))
offset(sol::Solution) = error("offset() not implemented for ", typeof(sol))

create_hdf5(sol::Solution, file) = error("create_hdf5() not implemented for ", typeof(sol))
CommonFunctions.write_to_hdf5(sol::Solution, h5::HDF5File, offset=0) = error("write_to_hdf5() not implemented for ", typeof(sol))

conv(sol::StochasticSolution) = error("conv() not implemented for ", typeof(sol))


"Create solution for ODE."
function Solution(equation::AbstractEquationODE, Δt, ntime::Int; kwargs...)
    SSolutionODE(equation, Δt, ntime; kwargs...)
end

"Create solution for partitioned ODE."
function Solution(equation::AbstractEquationPODE, Δt, ntime::Int; kwargs...)
    SSolutionPODE(equation, Δt, ntime; kwargs...)
end

"Create solution for DAE."
function Solution(equation::AbstractEquationDAE, Δt, ntime::Int; kwargs...)
    SSolutionDAE(equation, Δt, ntime; kwargs...)
end

"Create solution for partitioned DAE."
function Solution(equation::AbstractEquationPDAE, Δt, ntime::Int; kwargs...)
    SSolutionPDAE(equation, Δt, ntime; kwargs...)
end

"Create solution for SDE."
function Solution(equation::SDE, Δt, ntime::Int; kwargs...)
    SSolutionSDE(equation, Δt, ntime; kwargs...)
end

"Create solution for PSDE."
function Solution(equation::Union{PSDE,SPSDE}, Δt, ntime::Int; kwargs...)
    SSolutionPSDE(equation, Δt, ntime; kwargs...)
end

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, Δt, ntime::Int; kwargs...)
    error("No solution found for equation ", equation)
end

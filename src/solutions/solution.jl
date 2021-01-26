
using HDF5

abstract type Solution{dType, tType, N} end
abstract type DeterministicSolution{dType, tType, N} <: Solution{dType, tType, N} end

timesteps(sol::Solution) = error("time() not implemented for ", typeof(sol))
hdf5(sol::Solution)   = error("hdf5() not implemented for ", typeof(sol))
nsave(sol::Solution)  = error("nsave() not implemented for ", typeof(sol))
counter(sol::Solution) = error("counter() not implemented for ", typeof(sol))
offset(sol::Solution) = error("offset() not implemented for ", typeof(sol))
lastentry(sol::Solution) = error("lastentry() not implemented for ", typeof(sol))
Common.ntime(sol::Solution)  = error("ntime() not implemented for ", typeof(sol))

create_hdf5(sol::Solution, file) = error("create_hdf5() not implemented for ", typeof(sol))
Common.write_to_hdf5(sol::Solution, h5::HDF5.File, offset=0) = error("write_to_hdf5() not implemented for ", typeof(sol))

Common.nsamples(sol::DeterministicSolution) = sol.ni

Common.eachtimestep(sol::Solution) = 1:sol.nt*sol.nsave
Common.eachsample(sol::Solution) = 1:nsamples(sol)


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

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, Δt, ntime::Int; kwargs...)
    error("No solution found for equation ", equation)
end


"Create parallel solution for ODE."
function ParallelSolution(equation::AbstractEquationODE, Δt, ntime::Int; kwargs...)
    PSolutionODE(equation, Δt, ntime; kwargs...)
end

"Create parallel solution for partitioned ODE."
function ParallelSolution(equation::AbstractEquationPODE, Δt, ntime::Int; kwargs...)
    PSolutionPODE(equation, Δt, ntime; kwargs...)
end

"Create parallel solution for DAE."
function ParallelSolution(equation::AbstractEquationDAE, Δt, ntime::Int; kwargs...)
    PSolutionDAE(equation, Δt, ntime; kwargs...)
end

"Create parallel solution for partitioned DAE."
function ParallelSolution(equation::AbstractEquationPDAE, Δt, ntime::Int; kwargs...)
    PSolutionPDAE(equation, Δt, ntime; kwargs...)
end

"Print error for parallel solutions of equations not implemented, yet."
function ParallelSolution(equation::Equation, Δt, ntime::Int; kwargs...)
    error("No parallel solution found for equation ", equation)
end

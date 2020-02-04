
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
function Solution(equation::AbstractEquationODE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionODE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for partitioned ODE."
function Solution(equation::AbstractEquationPODE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionPODE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for DAE."
function Solution(equation::AbstractEquationDAE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionDAE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for partitioned DAE."
function Solution(equation::AbstractEquationPDAE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SSolutionPDAE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for SDE."
function Solution(equation::SDE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; K::Int=0, conv=:strong)
    SolutionSDE(equation, Δt, ntime, nsave, nwrite, K=K, conv=conv)
end

"Create solution for PSDE."
function Solution(equation::Union{PSDE,SPSDE}, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; K::Int=0, conv=:strong)
    SolutionPSDE(equation, Δt, ntime, nsave, nwrite, K=K, conv=conv)
end

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE)
    error("No solution found for equation ", equation)
end

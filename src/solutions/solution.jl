
using HDF5

abstract type Solution{dType, tType, N} end
abstract type DeterministicSolution{dType, tType, N} <: Solution{dType, tType, N} end
abstract type StochasticSolution{dType, tType, NQ, NW} <: Solution{dType, tType, NQ} end

hdf5(sol::Solution)   = error("hdf5() not implemented for ", typeof(sol))
time(sol::Solution)   = error("time() not implemented for ", typeof(sol))
ntime(sol::Solution)  = error("ntime() not implemented for ", typeof(sol))
nsave(sol::Solution)  = error("nsave() not implemented for ", typeof(sol))
offset(sol::Solution) = error("offset() not implemented for ", typeof(sol))

create_hdf5(sol::Solution, file) = error("create_hdf5() not implemented for ", typeof(sol))
CommonFunctions.write_to_hdf5(sol::Solution, h5::HDF5File, offset=0) = error("write_to_hdf5() not implemented for ", typeof(sol))


"Create solution for ODE and split ODE."
function Solution(equation::Union{ODE,SODE}, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionODE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for partitioned ODE."
function Solution(equation::PODE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionPODE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for variational ODE."
function Solution(equation::VODE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionPODE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for implicit ODE."
function Solution(equation::IODE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SSolutionPDAE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for DAE."
function Solution(equation::DAE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SolutionDAE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for partitioned DAE."
function Solution(equation::PDAE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SSolutionPDAE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for implicit DAE."
function Solution(equation::IDAE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE)
    SSolutionPDAE(equation, Δt, ntime, nsave, nwrite)
end

"Create solution for SDE."
function Solution(equation::SDE, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE; K::Int=0, conv::String="strong")
    SolutionSDE(equation, Δt, ntime, nsave, K=K, conv=conv)
end

"Create solution for PSDE."
function Solution(equation::Union{PSDE,SPSDE}, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE; K::Int=0, conv::String="strong")
    SolutionPSDE(equation, Δt, ntime, nsave, K=K, conv=conv)
end

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, Δt, ntime::Int, nsave::Int=DEFAULT_NSAVE)
    error("No solution found for equation ", equation)
end

"createHDF5: Creates or opens HDF5 file."
function createHDF5(sol::Solution, file::AbstractString; overwrite=true)
    if overwrite
        flag = "w"
        get_config(:verbosity) > 1 ? @info("Creating HDF5 file ", file) : nothing
        isfile(file) ? @warn("Overwriting existing HDF5 file.") : nothing
    else
        flag = "r+"
        get_config(:verbosity) > 1 ? @info("Opening HDF5 file ", file) : nothing
    end

    # create or open HDF5 file
    h5 = h5open(file, flag)

    return h5
end

"save_attributes: Saves common attributes of Solution to HDF5 file."
function save_attributes(sol::Solution)
    save_attributes(sol, hdf5(sol))
end

"save_attributes: Saves attributes of Deterministic Solutions to HDF5 file."
function save_attributes(sol::DeterministicSolution, h5::HDF5File)
    # save attributes
    attrs(h5)["ntime"] = ntime(sol)
    attrs(h5)["nsave"] = nsave(sol)
end

"save_attributes: Saves attributes of Stochastic Solutions to HDF5 file."
function save_attributes(sol::StochasticSolution, h5::HDF5File)
    attrs(h5)["nsave"] = solution.nsave
    attrs(h5)["conv"] = solution.conv
    attrs(h5)["nd"] = solution.nd
    attrs(h5)["nm"] = solution.nm
    attrs(h5)["ns"] = solution.ns
    attrs(h5)["ni"] = solution.ni
    attrs(h5)["K"] = solution.K
end

"write_to_hdf5: Wrapper for saving Solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::Solution)
    write_to_hdf5(solution, hdf5(solution), offset(solution))
end

"Creates HDF5 file, writes solution to file, and closes file."
function CommonFunctions.write_to_hdf5(solution::Solution, file::AbstractString)
    h5 = create_hdf5(solution, file)
    write_to_hdf5(solution, h5)
    close(h5)
end

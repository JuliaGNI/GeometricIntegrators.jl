
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
    attrs(h5)["ntime"] = ntime(sol)
    attrs(h5)["nsave"] = nsave(sol)
    attrs(h5)["conv"]  = string(conv(sol))
    attrs(h5)["nd"] = sol.nd
    attrs(h5)["nm"] = sol.nm
    attrs(h5)["ns"] = sol.ns
    attrs(h5)["K"]  = sol.K
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

# saving the time time series
function copy_timeteps_to_hdf5(sol::Solution, h5::HDF5File, j1, j2, n1, n2)
    h5["t"][j1:j2] = timesteps(sol)[n1:n2]
end

function copy_increments_to_hdf5(solution::StochasticSolution{DT,TT,NQ,2}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT,NQ}
    h5["ΔW"][:,j1:j2] = solution.W.ΔW[:,n1:n2]
    h5["ΔZ"][:,j1:j2] = solution.W.ΔZ[:,n1:n2]
end

function copy_increments_to_hdf5(solution::StochasticSolution{DT,TT,NQ,3}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT,NQ}
    h5["ΔW"][:,j1:j2,:] = solution.W.ΔW[:,n1:n2,:]
    h5["ΔZ"][:,j1:j2,:] = solution.W.ΔZ[:,n1:n2,:]
end

"""
Append solution to HDF5 file.
  offset - start writing q at the position offset+2
  offset2- start writing ΔW, ΔZ at the position offset2+1
"""
function CommonFunctions.write_to_hdf5(solution::StochasticSolution, h5::HDF5File=hdf5(solution), offset=offset(solution), offset2=offset)
    # set convenience variables and compute ranges
    d   = solution.nd
    m   = solution.nm
    s   = solution.ns

    j1  = offset+2
    j2  = offset+1+solution.nt
    jw1 = offset2+1
    jw2 = offset2+solution.ntime

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # saving the time time series
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    # copy data from solution to HDF5 dataset
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    if exists(h5, "ΔW") && exists(h5, "ΔZ")
        # copy the Wiener process increments from solution to HDF5 dataset
        copy_increments_to_hdf5(solution, h5, jw1, jw2, 1, solution.ntime)
    end

    return nothing
end

function determine_qdim(equation::Union{SDE,PSDE,SPSDE})
    ns = equation.ns
    ni = equation.ni

    @assert ns ≥ 1
    @assert ni ≥ 1

    if ns==ni==1
        NQ = 2
    elseif ns==1 || ni==1
        NQ = 3
    else
        @error("Both the number of sample paths and the number of initial conditions is larger than one!")
    end

    return NQ
end

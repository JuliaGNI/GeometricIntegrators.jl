
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
    # j1=offset(sol)+2, j2=offset(sol)+1+solution.nt, n1=1, n2=solution.nt
    set_dims!(h5["t"], (j2,))
    h5["t"][j1:j2] = timesteps(sol)[n1:n2]
end


function copy_solution_to_hdf5(solution::Union{SolutionODE{DT,TT,2},SolutionDAE{DT,TT,2},SolutionSDE{DT,TT,2}}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT}
    set_dims!(h5["q"], (size(h5["q"],1), j2))
    h5["q"][:, j1:j2] = solution.q[:, n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionODE{DT,TT,3},SolutionDAE{DT,TT,3},SolutionSDE{DT,TT,3}}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT}
    set_dims!(h5["q"], (size(h5["q"],1), j2, size(h5["q"],3)))
    h5["q"][:, j1:j2, :] = solution.q[:, n1:n2, :]
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{DT,TT,2},SolutionPDAE{DT,TT,2},SolutionPSDE{DT,TT,2}}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT}
    set_dims!(h5["q"], (size(h5["q"],1), j2))
    set_dims!(h5["p"], (size(h5["p"],1), j2))
    h5["q"][:, j1:j2] = solution.q[:, n1:n2]
    h5["p"][:, j1:j2] = solution.p[:, n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{DT,TT,3},SolutionPDAE{DT,TT,3},SolutionPSDE{DT,TT,3}}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT}
    set_dims!(h5["q"], (size(h5["q"],1), j2, size(h5["q"],3)))
    set_dims!(h5["p"], (size(h5["p"],1), j2, size(h5["p"],3)))
    h5["q"][:, j1:j2, :] = solution.q[:, n1:n2, :]
    h5["p"][:, j1:j2, :] = solution.p[:, n1:n2, :]
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{DT,TT,2},SolutionPDAE{DT,TT,2}}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT}
    set_dims!(h5["λ"], (size(h5["λ"],1), j2))
    h5["λ"][:, j1:j2] = solution.λ[:, n1:n2]
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{DT,TT,3},SolutionPDAE{DT,TT,3}}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT}
    set_dims!(h5["λ"], (size(h5["λ"],1), j2, size(h5["λ"],3)))
    h5["λ"][:, j1:j2, :] = solution.λ[:, n1:n2, :]
end

function copy_increments_to_hdf5(solution::StochasticSolution{DT,TT,NQ,2}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT,NQ}
    set_dims!(h5["ΔW"], (size(h5["ΔW"],1), j2))
    set_dims!(h5["ΔZ"], (size(h5["ΔZ"],1), j2))
    h5["ΔW"][:,j1:j2] = solution.W.ΔW[:,n1:n2]
    h5["ΔZ"][:,j1:j2] = solution.W.ΔZ[:,n1:n2]
end

function copy_increments_to_hdf5(solution::StochasticSolution{DT,TT,NQ,3}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT,NQ}
    set_dims!(h5["ΔW"], (size(h5["ΔW"],1), j2, size(h5["ΔW"],3)))
    set_dims!(h5["ΔZ"], (size(h5["ΔZ"],1), j2, size(h5["ΔZ"],3)))
    h5["ΔW"][:,j1:j2,:] = solution.W.ΔW[:,n1:n2,:]
    h5["ΔZ"][:,j1:j2,:] = solution.W.ΔZ[:,n1:n2,:]
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::DeterministicSolution, h5::HDF5File=hdf5(solution), offset=offset(solution))
    # set convenience variables and compute ranges
    j1 = offset+2
    j2 = offset+1+solution.nt

    # TODO # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::Union{SolutionDAE,SolutionPDAE}, h5::HDF5File=hdf5(solution), offset=offset(solution))
    # set convenience variables and compute ranges
    j1 = offset+2
    j2 = offset+1+solution.nt

    # TODO # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_multiplier_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    return nothing
end
"""
Append solution to HDF5 file.
  offset - start writing q at the position offset+2
  offset2- start writing ΔW, ΔZ at the position offset2+1
"""
function CommonFunctions.write_to_hdf5(solution::StochasticSolution, h5::HDF5File=hdf5(solution), offset=offset(solution), offset2=offset)
    # set convenience variables and compute ranges
    j1  = offset+2
    j2  = offset+1+solution.nt
    jw1 = offset2+1
    jw2 = offset2+solution.ntime

    # TODO # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    if exists(h5, "ΔW") && exists(h5, "ΔZ")
        # copy the Wiener process increments from solution to HDF5 dataset
        copy_increments_to_hdf5(solution, h5, jw1, jw2, 1, solution.ntime)
    end

    return nothing
end

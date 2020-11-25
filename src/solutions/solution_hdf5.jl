
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
function save_attributes(sol::DeterministicSolution, h5::HDF5.File)
    # save attributes
    attributes(h5)["ntime"] = ntime(sol)
    attributes(h5)["nsave"] = nsave(sol)
end

"save_attributes: Saves attributes of Stochastic Solutions to HDF5 file."
function save_attributes(sol::StochasticSolution, h5::HDF5.File)
    attributes(h5)["ntime"] = ntime(sol)
    attributes(h5)["nsave"] = nsave(sol)
    attributes(h5)["conv"]  = string(conv(sol))
    attributes(h5)["nd"] = sol.nd
    attributes(h5)["nm"] = sol.nm
    attributes(h5)["ns"] = sol.ns
    attributes(h5)["K"]  = sol.K
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



"Creates HDF5 file and initialises datasets for solution object."
function create_hdf5!(solution::Solution, file::AbstractString; kwargs...)
    solution.h5 = create_hdf5(solution, file; kwargs...)
end

"Creates HDF5 file and initialises datasets for deterministic solution object."
function create_hdf5(solution::DeterministicSolution, file::AbstractString)
    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution, h5)

    # create dataset
    init_timeteps_in_hdf5(solution, h5)
    init_solution_in_hdf5(solution, h5)

    return h5
end

"Creates HDF5 file and initialises datasets for stochastic solution object."
function create_hdf5(solution::StochasticSolution, file::AbstractString; save_W=true)
    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution, h5)

    # create dataset
    init_timeteps_in_hdf5(solution, h5)
    init_solution_in_hdf5(solution, h5)

    if save_W
        init_increments_in_hdf5(solution, h5)
    end

    return h5
end

function init_timeteps_in_hdf5(solution::Solution{DT,TT}, h5::HDF5.File) where {DT,TT}
    t = create_dataset(h5, "t", TT, ((solution.nt+1,), (-1,)), chunk=(1,))
    t[1] = timesteps(solution)[0]
end

function init_solution_in_hdf5(solution::SolutionODE{DT,TT,2}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    q[:, 1] = solution.q[:, 0]
end

function init_solution_in_hdf5(solution::SolutionDAE{DT,TT,2}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    λ = create_dataset(h5, "λ", DT, ((solution.nm, solution.nt+1), (solution.nm, -1)), chunk=(solution.nm,1))
    q[:, 1] = solution.q[:, 0]
    λ[:, 1] = solution.λ[:, 0]
end

function init_solution_in_hdf5(solution::SolutionSDE{DT,TT,2}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    q[:,1] = solution.q[:,0]
end

function init_solution_in_hdf5(solution::SolutionODE{DT,TT,3}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ni),(solution.nd, -1, solution.ni)), chunk=(solution.nd,1,1))
    q[:, 1, :] = solution.q[:, 0, :]
end

function init_solution_in_hdf5(solution::SolutionDAE{DT,TT,3}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ni),(solution.nd, -1, solution.ni)), chunk=(solution.nd,1,1))
    λ = create_dataset(h5, "λ", DT, ((solution.nm, solution.nt+1, solution.ni),(solution.nm, -1, solution.ni)), chunk=(solution.nm,1,1))
    q[:, 1, :] = solution.q[:, 0, :]
    λ[:, 1, :] = solution.λ[:, 0, :]
end

function init_solution_in_hdf5(solution::SolutionSDE{DT,TT,3}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ns),(solution.nd, -1, solution.ns)), chunk=(solution.nd,1,1))
    q[:,1,:] = solution.q[:,0,:]
end

function init_solution_in_hdf5(solution::SolutionPODE{DT,TT,2}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    q[:, 1] = solution.q[:, 0]
    p[:, 1] = solution.p[:, 0]
end

function init_solution_in_hdf5(solution::SolutionPDAE{DT,TT,2}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    λ = create_dataset(h5, "λ", DT, ((solution.nm, solution.nt+1), (solution.nm, -1)), chunk=(solution.nm,1))
    q[:, 1] = solution.q[:, 0]
    p[:, 1] = solution.p[:, 0]
    λ[:, 1] = solution.λ[:, 0]
end

function init_solution_in_hdf5(solution::SolutionPSDE{DT,TT,2}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    q[:,1] = solution.q[:,0]
    p[:, 1] = solution.p[:, 0]
end

function init_solution_in_hdf5(solution::SolutionPODE{DT,TT,3}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ni),(solution.nd, -1, solution.ni)), chunk=(solution.nd,1,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1, solution.ni),(solution.nd, -1, solution.ni)), chunk=(solution.nd,1,1))
    q[:, 1, :] = solution.q[:, 0, :]
    p[:, 1, :] = solution.p[:, 0, :]
end

function init_solution_in_hdf5(solution::SolutionPDAE{DT,TT,3}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ni),(solution.nd, -1, solution.ni)), chunk=(solution.nd,1,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1, solution.ni),(solution.nd, -1, solution.ni)), chunk=(solution.nd,1,1))
    λ = create_dataset(h5, "λ", DT, ((solution.nm, solution.nt+1, solution.ni),(solution.nm, -1, solution.ni)), chunk=(solution.nm,1,1))
    q[:, 1, :] = solution.q[:, 0, :]
    p[:, 1, :] = solution.p[:, 0, :]
    λ[:, 1, :] = solution.λ[:, 0, :]
end

function init_solution_in_hdf5(solution::SolutionPSDE{DT,TT,3}, h5::HDF5.File) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ns),(solution.nd, -1, solution.ns)), chunk=(solution.nd,1,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1, solution.ns),(solution.nd, -1, solution.ns)), chunk=(solution.nd,1,1))
    q[:,1,:] = solution.q[:,0,:]
    p[:, 1, :] = solution.p[:, 0, :]
end

function init_increments_in_hdf5(solution::StochasticSolution{DT,TT,NQ,2}, h5::HDF5.File) where {DT,TT,NQ}
    dW = create_dataset(h5, "ΔW", DT, ((solution.nm, solution.ntime),(solution.nm, -1)), chunk=(solution.nm,1))
    dZ = create_dataset(h5, "ΔZ", DT, ((solution.nm, solution.ntime),(solution.nm, -1)), chunk=(solution.nm,1))
end

function init_increments_in_hdf5(solution::StochasticSolution{DT,TT,NQ,3}, h5::HDF5.File) where {DT,TT,NQ}
    dW = create_dataset(h5, "ΔW", DT, ((solution.nm, solution.ntime, solution.ns),(solution.nm, -1, solution.ns)), chunk=(solution.nm,1,1))
    dZ = create_dataset(h5, "ΔZ", DT, ((solution.nm, solution.ntime, solution.ns),(solution.nm, -1, solution.ns)), chunk=(solution.nm,1,1))
end

function copy_timeteps_to_hdf5(sol::Solution, h5::HDF5.File, j1, j2, n1, n2)
    if size(h5["t"],1) < j2
        HDF5.set_dims!(h5["t"], (j2,))
    end
    h5["t"][j1:j2] = timesteps(sol)[n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionODE{DT,TT,2},SolutionDAE{DT,TT,2},SolutionSDE{DT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT}
    if size(h5["q"],2) < j2
        HDF5.set_dims!(h5["q"], (size(h5["q"],1), j2))
    end
    h5["q"][:, j1:j2] = solution.q[:, n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionODE{DT,TT,3},SolutionDAE{DT,TT,3},SolutionSDE{DT,TT,3}}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT}
    if size(h5["q"],2) < j2
        HDF5.set_dims!(h5["q"], (size(h5["q"],1), j2, size(h5["q"],3)))
    end
    h5["q"][:, j1:j2, :] = solution.q[:, n1:n2, :]
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{DT,TT,2},SolutionPDAE{DT,TT,2},SolutionPSDE{DT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT}
    if size(h5["q"],2) < j2
        HDF5.set_dims!(h5["q"], (size(h5["q"],1), j2))
    end
    if size(h5["p"],2) < j2
        HDF5.set_dims!(h5["p"], (size(h5["p"],1), j2))
    end
    h5["q"][:, j1:j2] = solution.q[:, n1:n2]
    h5["p"][:, j1:j2] = solution.p[:, n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{DT,TT,3},SolutionPDAE{DT,TT,3},SolutionPSDE{DT,TT,3}}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT}
    if size(h5["q"],2) < j2
        HDF5.set_dims!(h5["q"], (size(h5["q"],1), j2, size(h5["q"],3)))
    end
    if size(h5["p"],2) < j2
        HDF5.set_dims!(h5["p"], (size(h5["p"],1), j2, size(h5["p"],3)))
    end
    h5["q"][:, j1:j2, :] = solution.q[:, n1:n2, :]
    h5["p"][:, j1:j2, :] = solution.p[:, n1:n2, :]
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{DT,TT,2},SolutionPDAE{DT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT}
    if size(h5["λ"],2) < j2
        HDF5.set_dims!(h5["λ"], (size(h5["λ"],1), j2))
    end
    h5["λ"][:, j1:j2] = solution.λ[:, n1:n2]
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{DT,TT,3},SolutionPDAE{DT,TT,3}}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT}
    if size(h5["λ"],2) < j2
        HDF5.set_dims!(h5["λ"], (size(h5["λ"],1), j2, size(h5["λ"],3)))
    end
    h5["λ"][:, j1:j2, :] = solution.λ[:, n1:n2, :]
end

function copy_increments_to_hdf5(solution::StochasticSolution{DT,TT,NQ,2}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT,NQ}
    if size(h5["ΔW"],2) < j2
        HDF5.set_dims!(h5["ΔW"], (size(h5["ΔW"],1), j2))
    end
    if size(h5["ΔZ"],2) < j2
        HDF5.set_dims!(h5["ΔZ"], (size(h5["ΔZ"],1), j2))
    end
    h5["ΔW"][:,j1:j2] = solution.W.ΔW[:,n1:n2]
    h5["ΔZ"][:,j1:j2] = solution.W.ΔZ[:,n1:n2]
end

function copy_increments_to_hdf5(solution::StochasticSolution{DT,TT,NQ,3}, h5::HDF5.File, j1, j2, n1, n2) where {DT,TT,NQ}
    if size(h5["ΔW"],2) < j2
        HDF5.set_dims!(h5["ΔW"], (size(h5["ΔW"],1), j2, size(h5["ΔW"],3)))
    end
    if size(h5["ΔZ"],2) < j2
        HDF5.set_dims!(h5["ΔZ"], (size(h5["ΔZ"],1), j2, size(h5["ΔZ"],3)))
    end
    h5["ΔW"][:,j1:j2,:] = solution.W.ΔW[:,n1:n2,:]
    h5["ΔZ"][:,j1:j2,:] = solution.W.ΔZ[:,n1:n2,:]
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::DeterministicSolution, h5::HDF5.File=hdf5(solution), offset=offset(solution))
    # set convenience variables and compute ranges
    j1 = offset+2
    j2 = offset+1+solution.nt

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::Union{SolutionDAE,SolutionPDAE}, h5::HDF5.File=hdf5(solution), offset=offset(solution))
    # set convenience variables and compute ranges
    j1 = offset+2
    j2 = offset+1+solution.nt

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_multiplier_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    return nothing
end
"""
Append solution to HDF5 file.
  soffset - start writing the solution q at the position soffset+2
  woffset - start writing the increments ΔW, ΔZ at the position woffset+1
"""
function CommonFunctions.write_to_hdf5(solution::StochasticSolution, h5::HDF5.File=hdf5(solution), soffset=offset(solution), woffset=ioffset(solution))
    # set convenience variables and compute ranges
    js1 = soffset+2
    js2 = soffset+1+solution.nt
    jw1 = woffset+1
    jw2 = woffset+solution.nwrite

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, js1, js2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, js1, js2, 1, solution.nt)

    if haskey(h5, "ΔW") && haskey(h5, "ΔZ")
        # copy the Wiener process increments from solution to HDF5 dataset
        copy_increments_to_hdf5(solution, h5, jw1, jw2, 1, solution.nwrite)
    end

    return nothing
end

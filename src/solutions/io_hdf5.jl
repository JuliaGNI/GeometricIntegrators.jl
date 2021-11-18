"""

### Fields

* `h5`: HDF5 file for storage
* `nsave`: save data to disk after every nsave'th time step (default: ntime)
* `offset`: counter for file offset


"""
mutable struct SolutionHDF5 <: SolutionIO
    h5::HDF5.File
    nsave::Int
    offset::Int

    function SolutionHDF5(file::AbstractString; nsave::Int = 1, overwrite = false, readonly = false)
        if isfile(file)
            if readonly
                flag = "r"
                get_config(:verbosity) > 1 ? @info("Opening HDF5 file ", file) : nothing
            else
                if overwrite
                    flag = "w"
                    @warn("Overwriting existing HDF5 file ", file)
                else
                    flag = "r+"
                    get_config(:verbosity) > 1 ? @info("Opening HDF5 file ", file) : nothing
                end
            end
        else
            if readonly
                @error("Cannot open non-existing file in readonly mode.")
            else
                flag = "cw"
                get_config(:verbosity) > 1 ? @info("Creating HDF5 file ", file) : nothing
            end
        end

        # create or open HDF5 file
        h5 = h5open(file, flag)

        return new(h5, nsave, 0)
    end
end

function SolutionHDF5(file::AbstractString, solution::Solution; kwargs...)
    h5io = SolutionHDF5(file; kwargs...)
    initialize(h5io, solution)
    return h5io
end

hdf5(sio::SolutionHDF5) = sio.h5
GeometricBase.offset(sio::SolutionHDF5) = sio.offset
#attributes(sio::SolutionHDF5) = attributes(hdf5(sio))
Base.close(sio::SolutionHDF5) = close(hdf5(sio))


function initialize(h5::HDF5.File, solution::DeterministicSolution)
    # save attributes and common parameters
    save_attributes(h5, solution)

    # create dataset
    init_timeteps(h5, solution)
    init_solution(h5, solution)
end

"Save attributes and common parameters to HDF5 file and create data structures"
function initialize(sio::SolutionHDF5, solution::Solution)
    initialize(hdf5(sio), solution)
end


function save_attributes(h5::HDF5.File, solution::Solution)
    for (key, value) in pairs(parameters(solution))
        attributes(h5)[string(key)] = value
    end
end

"Saves attributes and parameters of a Solution to HDF5 file."
function save_attributes(sio::SolutionHDF5, solution::Solution)
    save_attributes(hdf5(sio), solution)
end


# function save(h5::HDF5.File, solution::Solution, offset)
#     # set convenience variables and compute ranges
#     j1 = offset + 2
#     j2 = offset + 1 + solution.nt

#     # copy data from solution to HDF5 dataset
#     save_timeteps(h5, solution, j1, j2, 1, solution.nt)
#     save_solution(h5, solution, j1, j2, 1, solution.nt)
#     save_multiplier(h5, solution, j1, j2, 1, solution.nt)

#     return nothing
# end

"Save Solution to HDF5 file."
function save(sio::SolutionHDF5, solution::Solution)
    # set convenience variables
    h5 = hdf5(sio)
    nt = solution.nt

    # compute ranges
    j1 = offset(sio) + 2
    j2 = offset(sio) + 1 + nt

    # copy data from solution to HDF5 dataset
    save_timeteps(h5, solution, j1, j2, 1, nt)
    save_solution(h5, solution, j1, j2, 1, nt)
    save_multiplier(h5, solution, j1, j2, 1, nt)

    sio.offset += nt

    return nothing
end


# function load(sio::SolutionHDF5)

# end



# "Creates HDF5 file, initialises datasets for deterministic solution object and runs user code."
# function create_hdf5(f::Function, file::AbstractString, solution::DeterministicSolution; kwargs...)
#     # create HDF5 file
#     h5 = _create_hdf5(file; kwargs...)

#     try
#         # save attributes and common parameters and create dataset
#         initialize_hdf5!(h5, solution)

#         # call user-provided code
#         f(h5)
#     finally
#         close(h5)
#     end
# end


# "Creates HDF5 file, writes solution to file, and closes file."
# function write_to_hdf5(file::AbstractString, solution::Solution)
#     create_hdf5(solution, file) do h5
#         write_to_hdf5(solution, h5)
#     end
# end



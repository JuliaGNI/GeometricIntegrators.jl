
using HDF5

# "createHDF5: Creates or opens HDF5 file."
# function createHDF5(file::AbstractString, overwrite=true)
#     if overwrite
#         flag = "w"
#     else
#         flag = "r+"
#     end
#
#     h5open(file, flag)
# end

"createHDF5: Creates HDF5 file and initialises datasets for solution object."
function createHDF5{T}(solution::SolutionODE{T}, file::AbstractString, ntime::Int=1)
    @assert ntime ≥ 1

    info("Creating HDF5 file ", file)
    # TODO Put warning if file exists.
    h5 = h5open(file, "w")

    # create dataset and copy initial conditions
    # ntime can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working.
    x = d_create(h5, "x", datatype(T), dataspace(solution.d, ntime))
    x[1:solution.d, 1] = solution[1:solution.d, 0]

    return h5
end

"writeSolutionToHDF5: Creates HDF5 file, writes solution to file, and closes file."
function writeSolutionToHDF5(solution::Solution, file::AbstractString)
    h5 = createHDF5(solution, file, solution.n+1)
    writeSolutionToHDF5(solution, h5)
    close(h5)
end

"writeSolutionToHDF5: Append solution to HDF5 file."
function writeSolutionToHDF5(solution::SolutionODE, h5::HDF5.HDF5File, offset=0)
    # aquire dataset from HDF5 file
    x = h5["x"]

    # set convenience variables and compute ranges
    d  = solution.d
    n  = solution.n
    j1 = offset+2
    j2 = offset+1+n

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    x[1:d, j1:j2] = solution[1:d, 1:n]

    return nothing
end

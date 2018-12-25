
using HDF5

abstract type Solution{dType, tType, N} end

abstract type StochasticSolution{dType, tType, NQ, NW} <: Solution{dType, tType, NQ} end

time(sol::Solution)  = error("time() not implemented for ", typeof(sol))
ntime(sol::Solution) = error("ntime() not implemented for ", typeof(sol))
nsave(sol::Solution) = error("nsave() not implemented for ", typeof(sol))


"Create solution for ODE and split ODE."
function Solution(equation::Union{ODE,SODE}, Δt, ntime::Int, nsave::Int=1)
    SolutionODE(equation, Δt, ntime, nsave)
end

"Create solution for partitioned ODE."
function Solution(equation::PODE, Δt, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, Δt, ntime, nsave)
end

"Create solution for variational ODE."
function Solution(equation::VODE, Δt, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, Δt, ntime, nsave)
end

"Create solution for implicit ODE."
function Solution(equation::IODE, Δt, ntime::Int, nsave::Int=1)
    SSolutionPDAE(equation, Δt, ntime, nsave)
end

"Create solution for DAE."
function Solution(equation::DAE, Δt, ntime::Int, nsave::Int=1)
    SolutionDAE(equation, Δt, ntime, nsave)
end

"Create solution for partitioned DAE."
function Solution(equation::PDAE, Δt, ntime::Int, nsave::Int=1)
    SSolutionPDAE(equation, Δt, ntime, nsave)
end

"Create solution for implicit DAE."
function Solution(equation::IDAE, Δt, ntime::Int, nsave::Int=1)
    SSolutionPDAE(equation, Δt, ntime, nsave)
end

"Create solution for SDE."
function Solution(equation::SDE, Δt, ntime::Int, nsave::Int=1; K::Int=0, conv::String="strong")
    SolutionSDE(equation, Δt, ntime, nsave, K=K, conv=conv)
end

"Create solution for PSDE."
function Solution(equation::Union{PSDE,SPSDE}, Δt, ntime::Int, nsave::Int=1; K::Int=0, conv::String="strong")
    SolutionPSDE(equation, Δt, ntime, nsave, K=K, conv=conv)
end

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, Δt, ntime::Int, nsave::Int=1)
    error("No solution found for equation ", equation)
end

"createHDF5: Creates or opens HDF5 file."
function createHDF5(sol::Solution, file::AbstractString, overwrite=true)
    if overwrite
        flag = "w"
        info("Creating HDF5 file ", file)
        isfile(file) ? warn("Overwriting existing HDF5 file.") : nothing
    else
        flag = "r+"
        info("Opening HDF5 file ", file)
    end

    # create or open HDF5 file
    h5 = h5open(file, flag)

    # save attributes
    attrs(h5)["ntime"] = ntime(sol)
    attrs(h5)["nsave"] = nsave(sol)

    # copy time
    write(h5, "t", time(sol))

    return h5
end


"createHDF5: Creates or opens HDF5 file."
# A version for StochasticSolution. It does not create attributes
# and does not write the time array t, like the version above does. Instead these
# are set in create_hdf5(), so that arrays larger than currently held in the solution
# structure can be created in the file. In the future it would be better to rewrite
# the function above, so that it is universal for all solution structures.
function createHDF5(sol::StochasticSolution, file::AbstractString, overwrite=true)
    if overwrite
        flag = "w"
        info("Creating HDF5 file ", file)
        isfile(file) ? warn("Overwriting existing HDF5 file.") : nothing
    else
        flag = "r+"
        info("Opening HDF5 file ", file)
    end

    # create or open HDF5 file
    h5 = h5open(file, flag)

    return h5
end


"Creates HDF5 file, writes solution to file, and closes file."
function writeSolutionToHDF5(solution::Solution, file::AbstractString)
    h5 = createHDF5(solution, file, solution.nt+1)
    writeSolutionToHDF5(solution, h5)
    close(h5)
end

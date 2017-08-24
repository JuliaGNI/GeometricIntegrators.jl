
using HDF5

abstract type Solution{dType, tType, N} end

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

"Creates HDF5 file, writes solution to file, and closes file."
function writeSolutionToHDF5(solution::Solution, file::AbstractString)
    h5 = createHDF5(solution, file, solution.nt+1)
    writeSolutionToHDF5(solution, h5)
    close(h5)
end


using HDF5

abstract Solution{dType, tType, N} <: AbstractArray{dType, N}

"Create solution for ODE."
function Solution(equation::ODE, Δt, ntime::Int, nsave::Int=1)
    SolutionODE(equation, Δt, ntime, nsave)
end

"Create solution for partitioned ODE."
function Solution(equation::PODE, Δt, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, Δt, ntime, nsave)
end

"Create solution for special ODE."
function Solution(equation::SODE, Δt, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, Δt, ntime, nsave)
end

"Create solution for DAE."
function Solution(equation::DAE, Δt, ntime::Int, nsave::Int=1)
    SolutionDAE(equation, ntime, nsave)
end

"Create solution for partitioned DAE."
function Solution(equation::PDAE, Δt, ntime::Int, nsave::Int=1)
    SolutionPDAE(equation, ntime, nsave)
end

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, Δt, ntime::Int, nsave::Int=1)
    error("No solution found for equation ", equation)
end

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

"Creates HDF5 file, writes solution to file, and closes file."
function writeSolutionToHDF5(solution::Solution, file::AbstractString)
    h5 = createHDF5(solution, file, solution.nt+1)
    writeSolutionToHDF5(solution, h5)
    close(h5)
end


Base.eltype{DT,TT,N}(s::Solution{DT,TT,N}) = DT
Base.ndims{DT,TT,N}(s::Solution{DT,TT,N}) = N
Base.size(s::Solution) = size(s.x)
# Base.length(s::Solution) = s.d * s.n
# Base.length(s::Solution) = length(s.x)
# Base.endof(s::Solution) = length(s)
Base.indices(s::Solution, d) = indices(s)[d]
Base.stride(s::Solution, d) = strides(s)[d]

# TODO Implement similar() and convert() to/from array functions.

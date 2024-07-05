
const default_extrapolation_stages = 5

abstract type Extrapolation end
# abstract type Extrapolation <: DeterministicMethod end


struct NoInitialGuess <: Extrapolation end


# """


# """
# function extrapolate! end



# function extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t₂, q₂, q̇₂, problem::AbstractProblemODE, extrap::Extrapolation; kwargs...)
#     solution = (
#         t = t₀, q = q₀, v = q̇₀,
#     )

#     history = (
#         (t = t₁, q = q₁, v = q̇₁),
#         (t = t₂, q = q₂, v = q̇₂),
#     )

#     solutionstep!(solution, history, problem, extrap; kwargs...)
# end

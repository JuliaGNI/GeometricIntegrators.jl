
# TODO (reactivate)

# n = 1
# T = Float64
#
# function F(x::Vector, b::Vector)
#     b[:] = x.^2
# end
#
# function J(x::Vector, A::Matrix)
#     A[:,:] = 2x
# end
#
#
# x = ones(T, n)
# nl = NewtonSolver(x, F)
# solve!(nl)
# # println(nl.status.i, ", ", nl.status.rₐ,", ",  nl.status.rᵣ,", ",  nl.status.rₛ)
# @test_approx_eq_eps(nl.x, zeros(T, n), 1E-7)
#
#
# x = ones(T, n)
# nl = NewtonSolver(x, F, J=J)
# solve!(nl)
# # println(nl.status.i, ", ", nl.status.rₐ,", ",  nl.status.rᵣ,", ",  nl.status.rₛ)
# @test_approx_eq_eps(nl.x, zeros(T, n), 1E-7)


n = 1
T = Float64

# J  = zeros(T, n, n)
# b  = zeros(T, n)
# lu = LUSolver(J, b)

function F(x::Vector, b::Vector)
    b[:] = x.^2
end

function J(x::Vector, A::Matrix)
    A[:,:] = 2x
end


x = ones(T, n)
nl = NewtonSolver(x, F)
solve!(nl)
# println(nl.i, ", ", nl.rₐ,", ",  nl.rᵣ,", ",  nl.rₛ)
@test_approx_eq_eps(nl.z, zeros(T, n), 1E-7)


x = ones(T, n)
nl = NewtonSolver(x, F, J=J)
solve!(nl)
# println(nl.i, ", ", nl.rₐ,", ",  nl.rᵣ,", ",  nl.rₛ)
@test_approx_eq_eps(nl.z, zeros(T, n), 1E-7)

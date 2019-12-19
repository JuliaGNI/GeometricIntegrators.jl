
using GeometricIntegrators.Solvers
using Test


n = 1
T = Float64

function F(x::Vector, b::Vector)
    b[:] = x.^2
end

function J(x::Vector, A::Matrix)
    A[:,:] = 2x
end


for Solver in (NewtonSolver, QuasiNewtonSolver)
    x = ones(T, n)
    nl = Solver(x, F)
    solve!(nl)
    # println(nl.status.i, ", ", nl.status.rₐ,", ",  nl.status.rᵣ,", ",  nl.status.rₛ)
    for x in nl.x
        @test x ≈ 0 atol=1E-7
    end

    x = ones(T, n)
    nl = Solver(x, F; J! = J)
    solve!(nl)
    # println(nl.status.i, ", ", nl.status.rₐ,", ",  nl.status.rᵣ,", ",  nl.status.rₛ)
    for x in nl.x
        @test x ≈ 0 atol=1E-7
    end
end

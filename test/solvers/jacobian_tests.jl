
using GeometricIntegrators.Config
using GeometricIntegrators.Solvers
using Test


n = 1
T = Float64

x = [T(π),]
j = reshape(2x, 1, 1)


function F(x::Vector, b::Vector)
    b[:] = x.^2
end

function J(x::Vector, A::Matrix)
    A[:,:] = 2x
end


set_config(:jacobian_autodiff, true)
JPAD = getJacobianParameters(nothing, F, T, n)

set_config(:jacobian_autodiff, false)
JPFD = getJacobianParameters(nothing, F, T, n)

set_config(:jacobian_autodiff, true)
JPUS = getJacobianParameters(J, F, T, n)

@test typeof(JPAD) <: JacobianParametersAD
@test typeof(JPFD) <: JacobianParametersFD
@test typeof(JPUS) <: JacobianParametersUser


function test_jac(j1, j2, atol)
    for i in eachindex(j1,j2)
        @test j1[i] ≈ j2[i] atol=atol
    end
end

jad = zero(j)
jfd = zero(j)
jus = zero(j)

computeJacobian(x, jad, JPAD)
computeJacobian(x, jfd, JPFD)
computeJacobian(x, jus, JPUS)

test_jac(jad, j, eps())
test_jac(jfd, j, 1E-7)
test_jac(jus, j, 0)

jad2 = zero(j)
jfd2 = zero(j)

computeJacobianAD(x, jad2, F)
computeJacobianFD(x, jfd2, F, get_config(:jacobian_fd_ϵ))

test_jac(jad, jad2, 0)
test_jac(jfd, jfd2, 0)

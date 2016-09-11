
abstract JacobianParameters{T}


immutable JacobianParametersUser{T} <: JacobianParameters{T}
    J::Function
end

immutable JacobianParametersAD{T} <: JacobianParameters{T}
    ϵ::T
    Fparams::NonlinearFunctionParameters{T}
end

immutable JacobianParametersFD{T} <: JacobianParameters{T}
    ϵ::T
    Fparams::NonlinearFunctionParameters{T}
    f1::Vector{T}
    f2::Vector{T}
    e::Vector{T}
    tx::Vector{T}
end


function computeJacobian{T}(x::Vector{T}, J::Matrix{T}, params::JacobianParametersUser{T})
    params.J(x, J)
end

function computeJacobian{T}(x::Vector{T}, J::Matrix{T}, params::JacobianParametersFD{T})
    @assert length(x) == size(J, 1) == size(J, 2)

    local ϵⱼ::T

    for j in 1:length(x)
        ϵⱼ = params.ϵ * x[j] + params.ϵ
        fill!(params.e, 0.)
        params.e[j] = 1.
        simd_waxpy!(params.tx, -ϵⱼ, params.e, x)
        function_stages!(params.tx, params.f1, params.Fparams)
        simd_waxpy!(params.tx, +ϵⱼ, params.e, x)
        function_stages!(params.tx, params.f2, params.Fparams)
        for i in 1:length(x)
            J[i,j] = (params.f2[i]-params.f1[i])/(2ϵⱼ)
        end
    end
end

# TODO (reactivate)
# function computeJacobianFD{T}(x::Vector{T}, J::Matrix{T}, F::Function, Fparams, ϵ::T)
#     local f1::Vector{T} = zeros(T, length(x))
#     local f2::Vector{T} = zeros(T, length(x))
#     local e::Vector{T}  = zeros(T, length(x))
#     local tx::Vector{T} = zeros(T, length(x))
#
#     computeJacobianFD(x, J, F, Fparams, ϵ, f1, f2, e, tx)
# end


function computeJacobian{T}(x::Vector{T}, J::Matrix{T}, Jparams::JacobianParametersAD{T})
    error("computeJacobian() with automatic differentiation not implemented, yet!")
end


function getJacobianParameters(J, Fparams, ϵ, T, n, autodiff)
    if J == nothing
        if autodiff
            Jparams = JacobianParametersAD{T}(ϵ, Fparams)
        else
            f1 = zeros(T, n)
            f2 = zeros(T, n)
            e  = zeros(T, n)
            tx = zeros(T, n)
            Jparams = JacobianParametersFD{T}(ϵ, Fparams, f1, f2, e, tx)
        end
    else
        @assert typeof(J) <: Function
        JacobianParametersUser{T}(J)
    end
    return Jparams
end

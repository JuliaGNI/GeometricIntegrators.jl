
using ForwardDiff


abstract JacobianParameters{T}


immutable JacobianParametersUser{T, FT <: Function} <: JacobianParameters{T}
    J::FT
end

immutable JacobianParametersAD{T, FT <: Function, FPT <: NonlinearFunctionParameters, N} <: JacobianParameters{T}
    ϵ::T
    f!::FT
    Fparams::FPT
    Jconfig::ForwardDiff.JacobianConfig{N}
    tx::Vector{T}
    ty::Vector{T}
end

function JacobianParametersAD{T, FT, FPT}(ϵ, f!::FT, Fparams::FPT, Jconfig, tx::Vector{T}, ty::Vector{T})
    @assert length(tx) == length(ty)
    JacobianParametersAD{T, FT, FPT, length(tx)}(ϵ, f!, Fparams, Jconfig, tx, ty)
end

immutable JacobianParametersFD{T, FPT <: NonlinearFunctionParameters} <: JacobianParameters{T}
    ϵ::T
    Fparams::FPT
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
        fill!(params.e, 0)
        params.e[j] = 1
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
    ForwardDiff.jacobian!(J, Jparams.f!, Jparams.ty, x, Jparams.Jconfig)
end


function getJacobianParameters(J, Fparams, ϵ, T, n, autodiff)
    if J == nothing
        if autodiff
            function function_stages_ad!(y, x)
                function_stages!(x, y, Fparams)
            end

            tx = zeros(T, n)
            ty = zeros(T, n)
            Jconfig = ForwardDiff.JacobianConfig(ty, tx)
            Jparams = JacobianParametersAD(ϵ, function_stages_ad!, Fparams, Jconfig, tx, ty)
        else
            f1 = zeros(T, n)
            f2 = zeros(T, n)
            e  = zeros(T, n)
            tx = zeros(T, n)
            Jparams = JacobianParametersFD{T, typeof(Fparams)}(ϵ, Fparams, f1, f2, e, tx)
        end
    else
        JacobianParametersUser{T, typeof(J)}(J)
    end
    return Jparams
end


using ForwardDiff


abstract type JacobianParameters{T} end


struct JacobianParametersUser{T, JT <: Function} <: JacobianParameters{T}
    J::JT
end

struct JacobianParametersAD{T, FT <: Function, N} <: JacobianParameters{T}
    F!::FT
    Jconfig::ForwardDiff.JacobianConfig{N}
    tx::Vector{T}
    ty::Vector{T}
end

function JacobianParametersAD(F!::FT, Jconfig, tx::Vector{T}, ty::Vector{T}) where {T, FT}
    @assert length(tx) == length(ty)
    JacobianParametersAD{T, FT, length(tx)}(F!, Jconfig, tx, ty)
end

function JacobianParametersAD(F!::FT, T, n::Int) where {FT <: Function}
    F!rev = (y,x) -> F!(x,y)
    tx = zeros(T, n)
    ty = zeros(T, n)
    Jconfig = ForwardDiff.JacobianConfig(nothing, ty, tx)
    Jparams = JacobianParametersAD(F!rev, Jconfig, tx, ty)
end


struct JacobianParametersFD{T, FT} <: JacobianParameters{T}
    ϵ::T
    F!::FT
    f1::Vector{T}
    f2::Vector{T}
    e::Vector{T}
    tx::Vector{T}
end

function JacobianParametersFD(F!::FT, ϵ, T, n::Int) where {FT <: Function}
    f1 = zeros(T, n)
    f2 = zeros(T, n)
    e  = zeros(T, n)
    tx = zeros(T, n)
    Jparams = JacobianParametersFD{T, typeof(F!)}(ϵ, F!, f1, f2, e, tx)
end


function computeJacobian(x::Vector{T}, J::Matrix{T}, params::JacobianParametersUser{T}) where {T}
    params.J(x, J)
end

function computeJacobian(x::Vector{T}, J::Matrix{T}, params::JacobianParametersFD{T}) where {T}
    @assert length(x) == size(J, 1) == size(J, 2)

    local ϵⱼ::T

    for j in 1:length(x)
        ϵⱼ = params.ϵ * x[j] + params.ϵ
        fill!(params.e, 0)
        params.e[j] = 1
        params.tx .= x .- ϵⱼ .* params.e
        params.F!(params.tx, params.f1)
        params.tx .= x .+ ϵⱼ .* params.e
        params.F!(params.tx, params.f2)
        for i in 1:length(x)
            J[i,j] = (params.f2[i]-params.f1[i])/(2ϵⱼ)
        end
    end
end

function computeJacobian(x::Vector{T}, J::Matrix{T}, Jparams::JacobianParametersAD{T}) where {T}
    ForwardDiff.jacobian!(J, Jparams.F!, Jparams.ty, x, Jparams.Jconfig)
end

function computeJacobianFD(x::Vector{T}, J::Matrix{T}, F!::FT, ϵ::T) where{T, FT <: Function}
    params = JacobianParametersFD(F!, ϵ, T, length(x))
    computeJacobian(x, J, params)
end

function computeJacobianAD(x::Vector{T}, J::Matrix{T}, F!::FT) where{T, FT <: Function}
    params = JacobianParametersAD(F!, T, length(x))
    computeJacobian(x, J, params)
end


function getJacobianParameters(J, F!, T, n)
    if J == nothing
        if get_config(:jacobian_autodiff)
            Jparams = JacobianParametersAD(F!, T, n)
        else
            Jparams = JacobianParametersFD(F!, get_config(:jacobian_fd_ϵ), T, n)
        end
    else
        JacobianParametersUser{T, typeof(J)}(J)
    end
    return Jparams
end



function computeJacobianFD{T}(x::Vector{T}, J::Matrix{T}, F::Function, ϵ::T)
    @assert length(x) == size(J, 1) == size(J, 2)

    local n = length(x)
    local ϵⱼ::T
    local e::Vector{T}
    local f1::Vector{T}
    local f2::Vector{T}

    f1 = zeros(T, n)
    f2 = zeros(T, n)

    for j in 1:length(x)
        ϵⱼ = ϵ * x[j] + ϵ
        e = zeros(n)
        e[j] = 1
        F(x - ϵⱼ*e, f1)
        F(x + ϵⱼ*e, f2)
        J[:,j] = (f2-f1)/(2ϵⱼ)
    end
end


function computeJacobianAD{T}(x::Vector{T}, J::Matrix{T}, F::Function, ϵ::T)
    error("computeJacobianAD() not implemented, yet!")
end


function getComputeJacobianFunction(J, F, ϵ, autodiff)
    if J == nothing
        if autodiff
            function computeJacobianADF(x, A)
                computeJacobianAD(x, A, F, ϵ)
            end
            J = computeJacobianADF
        else
            function computeJacobianFDF(x, A)
                computeJacobianFD(x, A, F, ϵ)
            end
            J = computeJacobianFDF
        end
    else
        @assert typeof(J) <: Function
    end
    return J
end

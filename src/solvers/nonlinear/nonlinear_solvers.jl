
abstract NonlinearSolver{T}

solve!(s::NonlinearSolver) = error("solve! not implemented for $(typeof(s))")


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


function residual_absolute{T}(x::Vector{T})
    local r::T = 0.
    for xᵢ in x
        r = max(r, xᵢ*xᵢ)
    end
    r
end

function residual_relative{T}(δx::Vector{T}, x::Vector{T})
    @assert length(x) == length(δx)
    local r::T = 0.
    for i in 1:length(x)
        r = max(r, abs(δx[i] / x[i]))
    end
    r
end

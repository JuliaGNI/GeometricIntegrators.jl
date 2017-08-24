

function compute_one_form(t, x, p, Θ)
    @assert size(x) == size(p)
    @assert length(t) == size(x,2) == size(p,2)

    tp = zeros(eltype(x), size(x,1))

    for k in 1:size(x,3)
        for j in 1:size(x,2)
            Θ(t[j], x[:,j,k], tp)
            p[:,j,k] .= tp
        end
    end
end


function compute_velocity(γ, γ̇)
    @assert size(γ) == size(γ̇)

    n = size(γ,2)
    k = 1im*2π*collect(0:div(n-1,2)+1)

    for l in 1:size(γ,1)
        γ̂ = rfft(γ[l,:])
        γ̇[l,:] .= irfft(k .* γ̂, n) / n
    end
end


function compute_correction(t, x, λ, σ, g)
    @assert size(x) == size(λ) == size(σ)
    @assert length(t) == size(x,2) == size(σ,2)

    tσ = zeros(eltype(x), size(x,1))

    for k in 1:size(x,3)
        for j in 1:size(x,2)
            g(t[j], x[:,j,k], λ[:,j,k], tσ)
            σ[:,j,k] .= tσ
        end
    end
end


function compute_loop_integral(p::Array{T,2}, v::Array{T,2}) where {T}
    local result = zero(T)
    local error  = zero(T)

    for i in 1:size(v,2)
        for l in 1:size(v,1)
            result, error = compensated_summation(p[l,i] * v[l,i], result, error)
        end
    end

    return result
end

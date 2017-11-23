
# general helper functions for Runge-Kutta methods.

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Matrix{T}, b::Vector{T}, Δt::T) where {T}
    @assert length(x) == length(xₑᵣᵣ)
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    for k in indices(ẋ, 1)
        for i in indices(ẋ, 2)
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[k,i], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::Union{Vector{T}, Vector{Double{T}}}, ẋ::Matrix{T}, b::Vector{T}, Δt::T) where {T}
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    local Δx::eltype(x)

    for k in indices(ẋ, 1)
        Δx = 0
        for i in indices(ẋ, 2)
            Δx += b[i] * ẋ[k,i]
        end
        x[k] += Δt * Δx
    end
end

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Matrix{T}, b::Vector{T}, b̂::Vector, Δt::T) where {T}
    update_solution!(x, xₑᵣᵣ, ẋ, b, Δt)
    update_solution!(x, xₑᵣᵣ, ẋ, b̂, Δt)
end

function update_solution!(x::Union{Vector{T}, Vector{Double{T}}}, ẋ::Matrix{T}, b::Vector{T}, b̂::Vector, Δt::T) where {T}
    update_solution!(x, ẋ, b, Δt)
    update_solution!(x, ẋ, b̂, Δt)
end


function cut_periodic_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, periodicity::Vector{T}) where {T}
    @assert length(x) == length(xₑᵣᵣ) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k] < 0
                (x[k], xₑᵣᵣ[k]) = compensated_summation(+periodicity[k], x[k], xₑᵣᵣ[k])
            end
            while x[k] ≥ periodicity[k]
                (x[k], xₑᵣᵣ[k]) = compensated_summation(-periodicity[k], x[k], xₑᵣᵣ[k])
            end
        end
    end
end

function cut_periodic_solution!(x::Vector{T}, periodicity::Vector{T}) where {T}
    @assert length(x) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k] < 0
                x[k] += periodicity[k]
            end
            while x[k] ≥ periodicity[k]
                x[k] -= periodicity[k]
            end
        end
    end
end

function cut_periodic_solution!(x::Vector{Double{T}}, periodicity::Vector{T}) where {T}
    @assert length(x) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k].hi < 0
                x[k] += periodicity[k]
            end
            while x[k].hi ≥ periodicity[k]
                x[k] -= periodicity[k]
            end
        end
    end
end

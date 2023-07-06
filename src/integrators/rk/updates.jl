
function update!(x::AbstractVector{T}, xₑᵣᵣ::AbstractVector{T}, ẋ::StageVector{T}, b::AbstractVector, Δt) where {T}
    @assert length(b) == length(ẋ)
    @assert length(x) == length(ẋ[1])
    @assert length(x) == length(xₑᵣᵣ)

    for k in eachindex(x,xₑᵣᵣ)
        for i in eachindex(ẋ)
            x[k], xₑᵣᵣ[k] = compensated_summation(x[k], Δt * b[i] * ẋ[i][k], xₑᵣᵣ[k])
        end
    end
end

function update!(x::AbstractVector{T}, xₑᵣᵣ::AbstractVector{T}, ẋ::StageVector{T}, b::AbstractVector, b̂::AbstractVector, Δt) where {T}
    update!(x, xₑᵣᵣ, ẋ, b, Δt)
    update!(x, xₑᵣᵣ, ẋ, b̂, Δt)
end

function update!(solstep::Union{SolutionStepODE,SolutionStepDAE}, V::StageVector, tableau::Tableau, Δt)
    update!(solstep.q, solstep.q̃, V, tableau.b, tableau.b̂, Δt)
end

function update!(solstep::Union{SolutionStepPODE,SolutionStepPDAE}, V::StageVector, F::StageVector, tableau::Tableau, Δt)
    update!(solstep.q, solstep.q̃, V, tableau.b, tableau.b̂, Δt)
    update!(solstep.p, solstep.p̃, F, tableau.b, tableau.b̂, Δt)
end

function update!(solstep::Union{SolutionStepPODE,SolutionStepPDAE}, V::StageVector, F::StageVector, tableau::PartitionedTableau, Δt)
    update!(solstep.q, solstep.q̃, V, tableau.q.b, tableau.q.b̂, Δt)
    update!(solstep.p, solstep.p̃, F, tableau.p.b, tableau.p.b̂, Δt)
end


function update_vector!(Δq::AbstractVector{T}, ẋ::StageVector{T}, b::AbstractVector, b̂::AbstractVector, Δt) where {T}
    @assert length(b) == length(ẋ)

    local δq₁::Base.TwicePrecision{T}
    local δq₂::Base.TwicePrecision{T}

    for k in eachindex(Δq)
        δq₁ = δq₂ = 0
        for i in eachindex(ẋ)
            δq₁ += Δt * b[i] * ẋ[i][k]
        end
        for i in eachindex(ẋ)
            δq₂ += Δt * b̂[i] * ẋ[i][k]
        end
        Δq[k] = (δq₁ + δq₂).hi
    end
end

function update_vector!(Δq::AbstractVector, V::StageVector, tableau::Tableau, Δt)
    update_vector!(Δq, V, tableau.b, tableau.b̂, Δt)
end

function update_vector!(Δq::AbstractVector, Δp::AbstractVector, V::StageVector, F::StageVector, tableau::Tableau, Δt)
    update_vector!(Δq, V, tableau.b, tableau.b̂, Δt)
    update_vector!(Δp, F, tableau.b, tableau.b̂, Δt)
end

function update_vector!(Δq::AbstractVector, Δp::AbstractVector, V::StageVector, F::StageVector, tableau::PartitionedTableau, Δt)
    update_vector!(Δq, V, tableau.q.b, tableau.q.b̂, Δt)
    update_vector!(Δp, F, tableau.p.b, tableau.p.b̂, Δt)
end


function update_multiplier!(λ::SolutionVector{T}, Λ::StageVector{T}, b::AbstractVector) where {T}
    for Λᵢ in Λ
        @assert length(λ) == length(Λᵢ)
    end
    local t::T
    @inbounds for i in eachindex(λ)
        t = zero(T)
        for j in eachindex(b,Λ)
            t += b[j] * Λ[j][i]
        end
        λ[i] = t
    end
    nothing
end


# function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Matrix{T}, b::AbstractVector, Δt) where {T}
#     @assert length(x) == length(xₑᵣᵣ)
#     @assert length(x) == size(ẋ, 1)
#     @assert length(b) == size(ẋ, 2)

#     for k in axes(ẋ, 1)
#         for i in axes(ẋ, 2)
#             x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[k,i], x[k], xₑᵣᵣ[k])
#         end
#     end
# end

# function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Vector{Vector{T}}, b::AbstractVector, Δt) where {T}
#     @assert length(b) == length(ẋ)
#     @assert length(x) == length(ẋ[1])
#     @assert length(x) == length(xₑᵣᵣ)

#     for i in eachindex(ẋ)
#         for k in eachindex(ẋ[i])
#             x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[i][k], x[k], xₑᵣᵣ[k])
#         end
#     end
# end

# function update_solution!(x::SolutionVector{T}, ẋ::Matrix{T}, b::AbstractVector, Δt) where {T}
#     @assert length(x) == size(ẋ, 1)
#     @assert length(b) == size(ẋ, 2)

#     local Δx::eltype(x)

#     for k in axes(ẋ, 1)
#         Δx = 0
#         for i in axes(ẋ, 2)
#             Δx += b[i] * ẋ[k,i]
#         end
#         x[k] += Δt * Δx
#     end
# end


# function update_solution!(x::SolutionVector{T}, ẋ::Vector{Vector{T}}, b::AbstractVector, Δt) where {T}
#     @assert length(b) == length(ẋ)
#     @assert length(x) == length(ẋ[1])

#     local Δx::T

#     for k in eachindex(x)
#         Δx = 0
#         for i in eachindex(ẋ)
#             Δx += b[i] * ẋ[i][k]
#         end
#         x[k] += Δt * Δx
#     end
# end

# function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::AbstractVector, b̂::AbstractVector, Δt) where {T}
#     update_solution!(x, xₑᵣᵣ, ẋ, b, Δt)
#     update_solution!(x, xₑᵣᵣ, ẋ, b̂, Δt)
# end

# function update_solution!(x::SolutionVector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::AbstractVector, b̂::AbstractVector, Δt) where {T}
#     update_solution!(x, ẋ, b, Δt)
#     update_solution!(x, ẋ, b̂, Δt)
# end

# function update_multiplier!(λ::SolutionVector{T}, Λ::Vector{Vector{T}}, b::AbstractVector) where {T}
#     for Λᵢ in Λ
#         @assert length(λ) == length(Λᵢ)
#     end
#     local t::T
#     @inbounds for i in eachindex(λ)
#         t = zero(T)
#         for j in eachindex(b,Λ)
#             t += b[j] * Λ[j][i]
#         end
#         λ[i] = t
#     end
#     nothing
# end

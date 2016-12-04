
type LUSolver{T} <: LinearSolver{T}
    n::Int
    A::Matrix{T}
    b::Vector{T}
    x::Vector{T}
    pivots::Vector{Int}
    info::Int

    LUSolver(n) = new(n, zeros(T, n, n), zeros(T, n), zeros(T, n), zeros(Int, n), 0)
end

function LUSolver{T}(A::Matrix{T})
    n = checksquare(A)
    lu = LUSolver{eltype(A)}(n)
    lu.A .= A
    lu
end

function LUSolver{T}(A::Matrix{T}, b::Vector{T})
    n = checksquare(A)
    @assert n == length(b)
    lu = LUSolver{eltype(A)}(n)
    lu.A .= A
    lu.b .= b
    lu
end


# function factorize!{T}(lu::LUSolver{T}, pivot=true)
#     @inbounds begin
#         for k = 1:lu.n
#             # find index max
#             kp = k
#             if pivot
#                 amax = real(zero(T))
#                 for i = k:lu.n
#                     absi = abs(lu.A[i,k])
#                     if absi > amax
#                         kp = i
#                         amax = absi
#                     end
#                 end
#             end
#             lu.pivots[k] = kp
#             if lu.A[kp,k] != 0
#                 if k != kp
#                     # Interchange
#                     for i = 1:lu.n
#                         tmp = lu.A[k,i]
#                         lu.A[k,i] = lu.A[kp,i]
#                         lu.A[kp,i] = tmp
#                     end
#                 end
#                 # Scale first column
#                 Akkinv = one(T)/lu.A[k,k]
#                 for i = k+1:lu.n
#                     lu.A[i,k] *= Akkinv
#                 end
#             elseif lu.info == 0
#                 lu.info = k
#             end
#             # Update the rest
#             for j = k+1:lu.n
#                 for i = k+1:lu.n
#                     lu.A[i,j] -= lu.A[i,k] * lu.A[k,j]
#                 end
#             end
#         end
#     end
#     nothing
# end


function factorize!{T}(lu::LUSolver{T})
    @inbounds for k = 1:lu.n
        lu.pivots[k] = k
        Ainv = one(T)/lu.A[k,k]
        @simd for i = k+1:lu.n
            lu.A[i,k] *= Ainv
            @simd for j = k+1:lu.n
                lu.A[i,j] -= lu.A[i,k] * lu.A[k,j]
            end
        end
    end
end


function solve!{T}(lu::LUSolver{T})
    for i = 1:lu.n
        lu.x[i] = lu.b[lu.pivots[i]]
    end
    for i = lu.n:-1:1
        for j = i+1:lu.n
            lu.x[i] -= lu.A[i,j] * lu.x[j]
        end
        lu.x[i] /= lu.A[i,i]
    end
    for i = 1:lu.n
        lu.b[i] = lu.x[i]
    end
    nothing
end

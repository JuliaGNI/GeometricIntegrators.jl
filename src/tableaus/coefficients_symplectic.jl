

function get_symplectic_conjugate_coefficients{T}(a::Matrix{T}, b::Vector{T}, a̅::Matrix{T})
    @assert size(a) == size(a̅)
    @assert length(b) == size(a,1) == size(a,2)

    for i in 1:size(a̅, 1)
        for j in 1:size(a̅, 2)
            a̅[i,j] = b[j] * ( one(T) - a[j,i] / b[i] )
        end
    end
    nothing
end


function get_symplectic_conjugate_coefficients{T}(coeff::CoefficientsRK{T})
    a̅ = zeros(coeff.a)
    get_symplectic_conjugate_coefficients(coeff.a, coeff.b, a̅)
    CoefficientsRK(coeff.name, coeff.s, a̅, coeff.b, coeff.c)
end

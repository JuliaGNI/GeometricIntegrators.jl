
abstract type AbstractCoefficients{T} end

Base.isequal(coeff1::AbstractCoefficients{T1}, coeff2::AbstractCoefficients{T2}) where {T1,T2} = (coeff1 == coeff2 && T1 == T2 && typeof(coeff1) == typeof(coeff2))

"""
Abstract basis

  T: data type
  N: number of nodes
"""
abstract type Basis{T<:AbstractFloat, N} end

nbasis(b::Basis{T,N}) where {T,N} = error("Not implemented!")
nnodes(b::Basis{T,N}) where {T,N} = error("Not implemented!")
degree(b::Basis{T,N}) where {T,N} = error("Not implemented!")
nodes(b::Basis{T,N})  where {T,N} = error("Not implemented!")

eval_basis(b::Basis, x, y)  = error("Not implemented!")
deriv_basis(b::Basis, x, y) = error("Not implemented!")
int_basis(b::Basis, x, y)   = error("Not implemented!")


"""
Abstract polynomial basis.
"""
abstract type PolynomialBasis{T,N} <: Basis{T,N} end

nbasis(b::PolynomialBasis{T,N}) where {T,N} = N
nnodes(b::PolynomialBasis{T,N}) where {T,N} = N
degree(b::PolynomialBasis{T,N}) where {T,N} = N-1

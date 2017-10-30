
abstract type Basis{T<:AbstractFloat, P} end

CommonFunctions.nbasis(b::Basis{T,P}) where {T,P} = error("Not implemented!")
CommonFunctions.nnodes(b::Basis{T,P}) where {T,P} = error("Not implemented!")
CommonFunctions.nodes(b::Basis{T,P})  where {T,P} = error("Not implemented!")
CommonFunctions.degree(b::Basis{T,P}) where {T,P} = error("Not implemented!")

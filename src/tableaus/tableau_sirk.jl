
"Holds the tableau of a singly implicit Runge-Kutta method."
immutable TableauSIRK{T} <: AbstractTableauIRK{T}
    # TODO
end

function TableauSIRK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauSIRK{T}(name, order, length(c), a, b, c)
end

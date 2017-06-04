
"Holds the tableau of a singly implicit Runge-Kutta method."
struct TableauSIRK{T} <: AbstractTableauIRK{T}
    # TODO
end

function TableauSIRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauSIRK{T}(name, order, length(c), a, b, c)
end

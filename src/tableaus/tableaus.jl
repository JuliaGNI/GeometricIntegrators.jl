
"Holds the information for the various methods' tableaus."
abstract type AbstractTableau{T} end


@define HeaderTableau begin
    name::Symbol
    o::Int
    s::Int
end

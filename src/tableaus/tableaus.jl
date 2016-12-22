
"Holds the information for the various methods' tableaus."
abstract AbstractTableau{T}


@define HeaderTableau begin
    name::Symbol
    o::Int
    s::Int
end

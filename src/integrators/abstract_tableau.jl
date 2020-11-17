
"Holds the information for the various methods' tableaus."
abstract type AbstractTableau{T} end

order(tab::T) where {T <: AbstractTableau} = tab.o
nstages(tab::T) where {T <: AbstractTableau} = tab.s

function check_symplecticity end
function symplecticity_conditions end

@define HeaderTableau begin
    name::Symbol
    o::Int
    s::Int
end

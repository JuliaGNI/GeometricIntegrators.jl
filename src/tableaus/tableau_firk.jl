
"Holds the tableau of a fully implicit Runge-Kutta method."
struct TableauFIRK{T} <: AbstractTableauIRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauFIRK{T}(q) where {T}
        if (q.s > 1 && istrilstrict(q.a)) || (q.s==1 && q.a[1,1] == 0)
            warn("Initializing TableauFIRK with explicit tableau ", q.name, ".\n",
                 "You might want to use TableauERK instead.")
        elseif q.s > 1 && istril(q.a)
            warn("Initializing TableauFIRK with diagonally implicit tableau ", q.name, ".\n",
                 "You might want to use TableauDIRK instead.")
        end

        new(q.name, q.o, q.s, q)
    end
end

function TableauFIRK(q::CoefficientsRK{T}) where {T}
    TableauFIRK{T}(q)
end

function TableauFIRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauFIRK{T}(CoefficientsRK(name, order, a, b, c))
end

# TODO function readTableauFIRKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a diagonally implicit Runge-Kutta method."
immutable TableauDIRK{T} <: AbstractTableauIRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauDIRK(q)
        @assert istril(q.a)
        @assert !(q.s==1 && q.a[1,1] ≠ 0)

        if q.s > 1 && istrilstrict(q.a)
            warn("Initializing TableauDIRK with explicit tableau ", q.name, ".\n",
                 "You might want to use TableauERK instead.")
        end

        new(q.name, q.o, q.s, q)
    end
end

function TableauDIRK{T}(q::CoefficientsRK{T})
    TableauDIRK{T}(q)
end

function TableauDIRK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauDIRK{T}(CoefficientsRK(name, order, a, b, c))
end

# TODO function readTableauDIRKFromFile(dir::AbstractString, name::AbstractString)

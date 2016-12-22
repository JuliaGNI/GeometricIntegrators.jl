
"Holds the tableau of an explicit Runge-Kutta method."
immutable TableauERK{T} <: AbstractTableauRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauERK(q)
        @assert q.c[1] == 0
        @assert istrilstrict(q.a)
        @assert !(q.s==1 && q.a[1,1] ≠ 0)
        new(q.name, q.o, q.s, q)
    end
end

function TableauERK{T}(q::CoefficientsRK{T})
    TableauERK{T}(q)
end

function TableauERK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauERK{T}(CoefficientsRK(name, order, a, b, c))
end

"Read explicit Runge-Kutta tableau from file."
function readTableauERKFromFile(dir::AbstractString, name::AbstractString)
    file = string(dir, "/", name, ".tsv")

#    run(`cat $file`)

    o, s, T = readTableauRKHeaderFromFile(file)

    # TODO Read data in original format (e.g., Rational).
    #      For this we need to save tableaus as jld or hdf5.
#    tab_array = readdlm(file, T)
    tab_array = readdlm(file)

    if s == 0
        s = size(tab_array, 1)
    end

    @assert s == size(tab_array, 1) == size(tab_array, 2)-2

    c = tab_array[1:s, 1]
    b = tab_array[1:s, 2]
    a = tab_array[1:s, 3:s+2]

    info("Reading explicit Runge-Kutta tableau ", name, " with ", s, " stages and order ", o, " from file\n      ", file)

    TableauERK(Symbol(name), o, a, b, c)
end

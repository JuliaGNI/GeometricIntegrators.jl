
"Holds the tableau of a Runge-Kutta method."
abstract type AbstractTableauRK{T <: Real} <: AbstractTableau{T} end

"Holds the tableau of an explicit Runge-Kutta method."
abstract type AbstractTableauERK{T} <: AbstractTableauRK{T} end

"Holds the tableau of an implicit Runge-Kutta method."
abstract type AbstractTableauIRK{T} <: AbstractTableauRK{T} end

"Holds the tableau of a partitioned Runge-Kutta method."
abstract type AbstractTableauPRK{T} <: AbstractTableauRK{T} end

Base.:(==)(tab1::T1, tab2::T2) where {T1 <: Union{AbstractTableauERK,AbstractTableauIRK},
                                      T2 <: Union{AbstractTableauERK,AbstractTableauIRK}} = (tab1.q == tab2.q && order(tab1) == order(tab2) && nstages(tab1) == nstages(tab2))

Base.:(==)(tab1::AbstractTableauPRK, tab2::AbstractTableauPRK) = (tab1.q == tab2.q && tab1.p == tab2.p && order(tab1) == order(tab2) && nstages(tab1) == nstages(tab2))

Base.isequal(tab1::AbstractTableauRK{DT1}, tab2::AbstractTableauRK{DT2}) where {DT1, DT2} = (tab1 == tab2 && DT1 == DT2 && typeof(tab1) == typeof(tab2))

"Reads and parses Tableau metadata from file."
function readTableauRKHeaderFromFile(file)
    f = open(file, "r")
    header = readline(f)
    close(f)

    if header[1] == '#'
        header = split(header[2:end])
    else
        header = ()
    end

    if length(header) ≥ 1
        O = Base.parse(Int, header[1])
    else
        O = 0
    end

    if length(header) ≥ 2
        S = Base.parse(Int, header[2])
    else
        S = 0
    end

    if length(header) ≥ 3
        T = Core.eval(Main, Meta.parse(header[3]))
    else
        T = Float64
    end

    return O, S, T
end

"Write Runge-Kutta tableau to file."
function writeTableauToFile(dir::AbstractString, tab::AbstractTableauRK{T}) where {T}
    # tab_array = zeros(T, S+1, S+1)
    # tab_array[1:S, 2:S+1] = tab.q.a
    # tab_array[S+1, 2:S+1] = tab.q.b
    # tab_array[1:S, 1] = tab.q.c
    # tab_array[S+1, 1] = tab.q.order

    tab_array = zeros(T, tab.q.s, tab.q.s+2)
    tab_array[1:tab.q.s, 1] = tab.q.c
    tab_array[1:tab.q.s, 2] = tab.q.b
    tab_array[1:tab.q.s, 3:tab.q.s+2] = tab.q.a

    header = string("# ", tab.q.o, " ", tab.q.s, " ", T, "\n")
    file   = string(dir, "/", tab.q.name, ".tsv")

    get_config(:verbosity) > 1 ? @info("Writing Runge-Kutta tableau $(tab.q.name) with $(tab.q.s) stages and order $(tab.q.o) to file\n$(file)") : nothing

    f = open(file, "w")
    write(f, header)
    writedlm(f, float(tab_array))
    close(f)
    # TODO Write data in original format (e.g., Rational).
end

# TODO function writeTableauToFile{Name, T}(dir::AbstractString, tab::AbstractTableauPRK{Name, T})

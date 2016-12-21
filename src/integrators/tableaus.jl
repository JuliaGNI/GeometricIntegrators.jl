
"Holds the information for the various methods' tableaus."
abstract Tableau{T}

"Holds the tableau of a Runge-Kutta method."
abstract TableauRK{T} <: Tableau{T}

"Holds the tableau of an implicit Runge-Kutta method."
abstract TableauIRK{T} <: TableauRK{T}

"Holds the tableau of a partitioned Runge-Kutta method."
abstract TableauPRK{T} <: TableauRK{T}

@define HeaderTableauRK begin
    name::Symbol
    o::Int
    s::Int
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}
end

@define HeaderTableauPRK begin
    name::Symbol
    o::Int
    s::Int
    a_q::Matrix{T}
    a_p::Matrix{T}
    b_q::Vector{T}
    b_p::Vector{T}
    c_q::Vector{T}
    c_p::Vector{T}
end

Base.hash(tab::TableauRK, h::UInt) = hash(tab.o, hash(tab.a, hash(tab.b, hash(tab.c, hash(:TableauRK, h)))))
Base.:(==){T1, T2}(tab1::TableauRK{T1}, tab2::TableauRK{T2}) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && tab1.a == tab2.a
                                                             && tab1.b == tab2.b
                                                             && tab1.c == tab2.c)

Base.isequal{T1, T2}(tab1::TableauRK{T1}, tab2::TableauRK{T2}) = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

"Print Runge-Kutta tableau to standard output."
function showTableau{T}(tab::TableauRK{T})
    println("Runge-Kutta Method ", tab.name, "with ", tab.s, " stages and order ", tab.o)
    println("  a = ", tab.a)
    println("  b = ", tab.b)
    println("  c = ", tab.c)
end

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
        O = parse(Int, header[1])
    else
        O = 0
    end

    if length(header) ≥ 2
        S = parse(Int, header[2])
    else
        S = 0
    end

    if length(header) ≥ 3
        T = eval(parse(header[3]))
    else
        T = Float64
    end

    return O, S, T
end

"Write Runge-Kutta tableau to file."
function writeTableauToFile{T}(dir::AbstractString, tab::TableauRK{T})
    # tab_array = zeros(T, S+1, S+1)
    # tab_array[1:S, 2:S+1] = tab.a
    # tab_array[S+1, 2:S+1] = tab.b
    # tab_array[1:S, 1] = tab.c
    # tab_array[S+1, 1] = tab.order

    tab_array = zeros(T, tab.s, tab.s+2)
    tab_array[1:tab.s, 1] = tab.c
    tab_array[1:tab.s, 2] = tab.b
    tab_array[1:tab.s, 3:tab.s+2] = tab.a

    header = string("# ", tab.o, " ", tab.s, " ", T, "\n")
    file   = string(dir, "/", tab.name, ".tsv")

    info("Writing Runge-Kutta tableau ", tab.name, " with ", tab.s, " stages and order ", tab.o, " to file\n      ", file)

    f = open(file, "w")
    write(f, header)
    writedlm(f, float(tab_array))
    close(f)
    # TODO Write data in original format (e.g., Rational).
end

# TODO function writeTableauToFile{Name, T}(dir::AbstractString, tab::TableauPRK{Name, T})
# TODO function writeTableauToFile{Name, T}(dir::AbstractString, tab::TableauSPARK{Name, T})
# TODO function writeTableauToFile{Name, T}(dir::AbstractString, tab::TableauGLM{Name, T})


"Holds the tableau of an explicit Runge-Kutta method."
immutable TableauERK{T} <: TableauRK{T}
    @HeaderTableauRK

    function TableauERK(name,o,s,a,b,c)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert c[1] == 0
        @assert istrilstrict(a)
        @assert !(s==1 && a[1,1] ≠ 0)
        new(name,o,s,a,b,c)
    end
end

function TableauERK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauERK{T}(name, order, length(c), a, b, c)
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


"Holds the tableau of a diagonally implicit Runge-Kutta method."
immutable TableauDIRK{T} <: TableauIRK{T}
    @HeaderTableauRK

    function TableauDIRK(name,o,s,a,b,c)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert istril(a)
        @assert !(s==1 && a[1,1] ≠ 0)

        if s > 1 && istrilstrict(a)
            warn("Initializing TableauDIRK with explicit tableau ", name, ".\n",
                 "You might want to use TableauERK instead.")
        end

        new(name,o,s,a,b,c)
    end
end

function TableauDIRK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauDIRK{T}(name, order, length(c), a, b, c)
end

# TODO function readTableauDIRKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a fully implicit Runge-Kutta method."
immutable TableauFIRK{T} <: TableauIRK{T}
    @HeaderTableauRK

    function TableauFIRK(name,o,s,a,b,c)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)

        if (s > 1 && istrilstrict(a)) || (s==1 && a[1,1] == 0)
            warn("Initializing TableauFIRK with explicit tableau ", name, ".\n",
                 "You might want to use TableauERK instead.")
        elseif s > 1 && istril(a)
            warn("Initializing TableauFIRK with diagonally implicit tableau ", name, ".\n",
                 "You might want to use TableauDIRK instead.")
        end

        new(name,o,s,a,b,c)
    end
end

function TableauFIRK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauFIRK{T}(name, order, length(c), a, b, c)
end

# TODO function readTableauFIRKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a singly implicit Runge-Kutta method."
immutable TableauSIRK{T} <: TableauIRK{T}
    # TODO
end

function TableauSIRK{T}(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauSIRK{T}(name, order, length(c), a, b, c)
end


"Holds the tableau of an explicit partitioned Runge-Kutta method."
# TODO Need explicit and implicit version?
immutable TableauEPRK{T} <: Tableau{T}
    @HeaderTableauPRK

    function TableauEPRK(name, o, s, a_q, a_p, b_q, b_p, c_q, c_p)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        # TODO check that both tableaus are lower triangular and that only one element
        #      a_q[i,i] or a_p[i,i] is non-zero for all i.
        new(name, o, s, a_q, a_p, b_q, b_p, c_q, c_p)
    end
end

function TableauEPRK{T}(name::Symbol, order::Int,
                        a_q::Matrix{T}, a_p::Matrix{T},
                        b_q::Vector{T}, b_p::Vector{T},
                        c_q::Vector{T}, c_p::Vector{T})
    @assert length(c_q) == length(c_p)
    TableauEPRK{T}(name, order, length(c_q), a_q, a_p, b_q, b_p, c_q, c_p)
end

function TableauEPRK{T}(name::Symbol, order::Int, tab::TableauERK{T})
    TableauEPRK{T}(name, order, tab.s, tab.a, tab.a, tab.b, tab.b, tab.c, tab.c)
end

function TableauEPRK{T}(name::Symbol, order::Int, tab_q::TableauRK{T}, tab_p::TableauRK{T})
    @assert tab_q.s == tab_p.s
    TableauEPRK{T}(name, order, tab_q.s, tab_q.a, tab_p.a, tab_q.b, tab_p.b, tab_q.c, tab_p.c)
end

# TODO function readTableauPRKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a implicit partitioned Runge-Kutta method."
immutable TableauIPRK{T} <: Tableau{T}
    @HeaderTableauPRK

    function TableauIPRK(name, o, s, a_q, a_p, b_q, b_p, c_q, c_p)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)

        new(name, o, s, a_q, a_p, b_q, b_p, c_q, c_p)
    end
end

function TableauIPRK{T}(name::Symbol, order::Int,
                        a_q::Matrix{T}, a_p::Matrix{T},
                        b_q::Vector{T}, b_p::Vector{T},
                        c_q::Vector{T}, c_p::Vector{T})
    @assert length(c_q) == length(c_p)
    TableauIPRK{T}(name, order, length(c_q), a_q, a_p, b_q, b_p, c_q, c_p)
end

function TableauIPRK{T}(name::Symbol, order::Int, tab_q::TableauRK{T}, tab_p::TableauRK{T})
    @assert tab_q.s == tab_p.s
    TableauIPRK{T}(name, order, tab_q.s, tab_q.a, tab_p.a, tab_q.b, tab_p.b, tab_q.c, tab_p.c)
end

# TODO function readTableauIPRKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a variational partitioned Runge-Kutta method."
immutable TableauVPRK{T} <: Tableau{T}
    @HeaderTableauPRK

    d::Vector{T}

    function TableauVPRK(name, o, s, a_q, a_p, b_q, b_p, c_q, c_p, d)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages must be > 0"
        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        @assert s==length(d)

        new(name, o, s, a_q, a_p, b_q, b_p, c_q, c_p, d)
    end
end

function TableauVPRK{T}(name::Symbol, order::Int,
                        a_q::Matrix{T}, a_p::Matrix{T},
                        b_q::Vector{T}, b_p::Vector{T},
                        c_q::Vector{T}, c_p::Vector{T},
                        d::Vector{T})
    @assert length(c_q) == length(c_p)
    TableauVPRK{T}(name, order, length(c_q), a_q, a_p, b_q, b_p, c_q, c_p, d)
end

function TableauVPRK{T}(name::Symbol, order::Int, tab_q::TableauRK{T}, tab_p::TableauRK{T}, d::Vector{T})
    @assert tab_q.s == tab_p.s
    TableauVPRK{T}(name, order, tab_q.s, tab_q.a, tab_p.a, tab_q.b, tab_p.b, tab_q.c, tab_p.c, d)
end

# TODO function readTableauVPRKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a spezialized additive Runge-Kutta method."
immutable TableauSARK{T} <: Tableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    a_q::Matrix{T}
    α_q::Matrix{T}

    a_q̃::Matrix{T}
    α_q̃::Matrix{T}

    b_q::Vector{T}
    β_q::Vector{T}

    c_q::Vector{T}
    c_λ::Vector{T}

    ω_q::Matrix{T}
    ω_λ::Matrix{T}

    function TableauSARK(name, o, s, r,
                         a_q, α_q, a_q̃, α_q̃,
                         b_q, β_q, c_q, c_λ,
                         ω_q, ω_λ)
        # TODO Make ω_q, ω_λ optional arguments.
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(α_q,1)==length(β_q)
        @assert r==size(α_q,2)==size(α_p,2)
        @assert r==length(c_λ)
        @assert r==size(a_qᵠ,1)==size(α_qᵠ,1)==size(α_qᵠ,2)
        @assert s==size(a_qᵠ,2)
        # TODO Add assertions on ω_q, ω_λ to be (S-1)x(S) or (R-1)x(R) if set.
        new(name, o, s, r, a_q, α_q, a_q̃, α_q̃, b_q, β_q, c_q, c_λ, ω_q, ω_λ)
    end
end

# TODO Add external constructor for TableauSARK.

# TODO function readTableauSARKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a spezialized partitioned additive Runge-Kutta method."
immutable TableauSPARK{T} <: Tableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    a_q::Matrix{T}
    a_p::Matrix{T}
    α_q::Matrix{T}
    α_p::Matrix{T}

    a_q̃::Matrix{T}
    a_p̃::Matrix{T}
    α_q̃::Matrix{T}
    α_p̃::Matrix{T}

    b_q::Vector{T}
    b_p::Vector{T}
    β_q::Vector{T}
    β_p::Vector{T}

    c_q::Vector{T}
    c_p::Vector{T}
    c_λ::Vector{T}

    ω_q::Matrix{T}
    ω_p::Matrix{T}
    ω_λ::Matrix{T}

    function TableauSPARK(name, o, s, r,
                          a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                          b_q, b_p, β_q, β_p,
                          c_q, c_p, c_λ,
                          ω_q, ω_p, ω_λ)
        # TODO Make ω_q, ω_p, ω_λ optional arguments.
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        @assert s==size(α_q,1)==size(α_p,1)==length(β_q)==length(β_p)
        @assert r==size(α_q,2)==size(α_p,2)
        @assert r==length(c_λ)
        @assert r==size(a_qᵠ,1)==size(a_pᵠ,1)
        @assert r==size(α_qᵠ,1)==size(α_qᵠ,2)
        @assert r==size(α_pᵠ,1)==size(α_pᵠ,2)
        @assert s==size(a_qᵠ,2)==size(a_pᵠ,2)
        # TODO Add assertions on ω_q, ω_p, ω_λ to be (S-1)x(S) or (R-1)x(R) if set.
        new(name, o, s, r,
            a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
            b_q, b_p, β_q, β_p,
            c_q, c_p, c_λ,
            ω_q, ω_p, ω_λ)
    end
end

# TODO Add external constructor for TableauSPARK.

# TODO function readTableauSPARKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of a general linear method."
immutable TableauGLM{T} <: Tableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    a::Matrix{T}
    b::Matrix{T}
    u::Matrix{T}
    v::Matrix{T}
    c::Vector{T}

    function TableauGLM(name, o, s, r, a, b, u, v, c)
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(s, Integer)
        @assert isa(r, Integer)
        @assert isa(o, Integer)
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert s==length(c)
        @assert s==size(a,1)==size(a,2)==length(u,1)==length(b,2)
        @assert r==size(v,1)==size(v,2)==length(u,2)==length(b,1)
        new(name, o, s, r, a, b, u, v, c)
    end
end

# TODO Add external constructor for TableauGLM.

# TODO function readTableauGLMFromFile(dir::AbstractString, name::AbstractString)


# TODO Add TableauSGALM. # Super General Additive Linear Method
# TODO Add TableauTSRK.  # Two Step Runge Kutta Method
# TODO Add TableauSplitting.

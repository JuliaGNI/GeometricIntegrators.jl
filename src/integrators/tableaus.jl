
"Tableau: Holds the information for the various methods' tableaus."
abstract Tableau{Name, T}

"TableauRK: Holds the tableau of a Runge-Kutta method."
abstract TableauRK{Name, S, T} <: Tableau{Name, T}

Base.hash(tab::TableauRK, h::UInt) = hash(tab.order, hash(tab.a, hash(tab.b, hash(tab.c, hash(:TableauRK, h)))))
Base.(:(==)){Name1, Name2, S1, S2, T1, T2}(tab1::TableauRK{Name1, S1, T1}, tab2::TableauRK{Name2, S2, T2}) = (tab1.order == tab2.order
                                               && tab1.a == tab2.a
                                               && tab1.b == tab2.b
                                               && tab1.c == tab2.c)

Base.isequal{Name1, Name2, S1, S2, T1, T2}(tab1::TableauRK{Name1, S1, T1}, tab2::TableauRK{Name2, S2, T2}) = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

function showTableau{Name, S, T}(tab::TableauRK{Name, S, T})
    println("Explicit Runge-Kutta Method ", Name, "with ", S, " stages and order ", tab.order)
    println("  a = ", tab.a)
    println("  b = ", tab.b)
    println("  c = ", tab.c)
end

"readTableauHeaderFromFile: Reads Tableau metadata from file."
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

function writeTableauToFile{Name, S, T}(dir::AbstractString, tab::TableauRK{Name, S, T})
    # tab_array = zeros(T, S+1, S+1)
    # tab_array[1:S, 2:S+1] = tab.a
    # tab_array[S+1, 2:S+1] = tab.b
    # tab_array[1:S, 1] = tab.c
    # tab_array[S+1, 1] = tab.order

    tab_array = zeros(T, S, S+2)
    tab_array[1:S, 1] = tab.c
    tab_array[1:S, 2] = tab.b
    tab_array[1:S, 3:S+2] = tab.a

    header = string("# ", tab.order, " ", S, " ", T, "\n")
    file   = string(dir, "/", Name, ".tsv")

    println("  Writing Runge-Kutta tableau ", Name, " with ", S, " stages and order ", tab.order, " to file")
    println("  ", file, ".")

    f = open(file, "w")
    write(f, header)
    writedlm(f, float(tab_array))
    close(f)
    # TODO Write data in original format (e.g., Rational).
end

# TODO function writeTableauToFile{Name, S, T}(dir::AbstractString, tab::TableauPRK{Name, S, T})
# TODO function writeTableauToFile{Name, S, T}(dir::AbstractString, tab::TableauSPARK{Name, S, R, T})
# TODO function writeTableauToFile{Name, S, T}(dir::AbstractString, tab::TableauGLM{Name, S, R, T})


"TableauERK: Holds the tableau of an explicit Runge-Kutta method."
immutable TableauERK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function TableauERK(order,a,b,c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S > 0
        @assert S==size(a,1)==size(a,2)==length(b)==length(c)
        @assert c[1] == 0
        @assert istrilstrict(a)
        new(order,a,b,c)
    end
end

function TableauERK{T}(name::Symbol, order::Integer,
                       a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauERK{name, length(c), T}(order, a, b, c)
end

function readTableauERKFromFile(dir::AbstractString, name::AbstractString)
    file = string(dir, "/", name, ".tsv")

#    run(`cat $file`)

    order, S, T = readTableauRKHeaderFromFile(file)

    # TODO Read data in original format (e.g., Rational).
#    tab_array = readdlm(file, T)
    tab_array = readdlm(file)

    if S == 0
        S = size(tab_array, 1)
    end

    @assert S == size(tab_array, 1) == size(tab_array, 2)-2

    c = tab_array[1:S, 1]
    b = tab_array[1:S, 2]
    a = tab_array[1:S, 3:S+2]

    println("  Reading explicit Runge-Kutta tableau ", name, " with ", S, " stages and order ", order, " from file")
    println("  ", file)

    TableauERK(symbol(name), order, a, b, c)
end


"TableauIRK: Holds the tableau of a linearly implicit Runge-Kutta method."
immutable TableauIRK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function TableauIRK(order,a,b,c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S > 0
        @assert S==size(a,1)==size(a,2)==length(b)==length(c)
        @assert istril(a)

        if S == 1 # catch the case of a 1x1 matrix a
        elseif istrilstrict(a)
            warn("Initializing TableauIRK with explicit tableau.")
            info("   You might want to use TableauERK instead.")
        end

        new(order,a,b,c)
    end
end

function TableauIRK{T}(name::Symbol, order::Integer,
                       a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauIRK{name, length(c), T}(order, a, b, c)
end

# TODO function readTableauIRKFromFile(dir::AbstractString, name::AbstractString)


"TableauNLIRK: Holds the tableau of a nonlinearly implicit Runge-Kutta method."
immutable TableauNLIRK{Name, S, T} <: TableauRK{Name, S, T}
    order::Integer
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function TableauNLIRK(order,a,b,c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S > 0
        @assert S==size(a,1)==size(a,2)==length(b)==length(c)

        if S == 1 # catch the case of a 1x1 matrix a
        elseif istrilstrict(a)
            warn("Initializing TableauNLIRK with explicit tableau.")
            info("   You might want to use TableauERK instead.")
        elseif istril(a)
            warn("Initializing TableauNLIRK with linearly implicit tableau.")
            info("   You might want to use TableauIRK instead.")
        end

        new(order,a,b,c)
    end
end

function TableauNLIRK{T}(name::Symbol, order::Integer,
                         a::Matrix{T}, b::Vector{T}, c::Vector{T})
    TableauNLIRK{name, length(c), T}(order, a, b, c)
end

# TODO function readTableauNLIRKFromFile(dir::AbstractString, name::AbstractString)


"TableauPRK: Holds the tableau of a partitioned Runge-Kutta method."
# TODO Need explicit and implicit version?
immutable TableauPRK{Name, S, T} <: Tableau{Name, T}
    order::Integer
    a_q::Matrix{T}
    a_p::Matrix{T}
    b_q::Vector{T}
    b_p::Vector{T}
    c_q::Vector{T}
    c_p::Vector{T}

    function TableauPRK(order, a_q, a_p, b_q, b_p, c_q, c_p)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(order, Integer)
        @assert S > 0
        @assert S==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert S==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        new(order, a_q, a_p, b_q, b_p, c_q, c_p)
    end
end

function TableauPRK{T}(name::Symbol, order::Integer,
                       a_q::Matrix{T}, a_p::Matrix{T},
                       b_q::Vector{T}, b_p::Vector{T},
                       c_q::Vector{T}, c_p::Vector{T})
    @assert length(c_q)==length(c_p)
    TableauPRK{name, length(c_q), T}(order, a_q, a_p, b_q, b_p, c_q, c_p)
end

# TODO function readTableauPRKFromFile(dir::AbstractString, name::AbstractString)


"TableauSPARK: Holds the tableau of a spezialized partitioned additive
 Runge-Kutta method."
immutable TableauSPARK{Name, S, R, T} <: Tableau{Name, T}
    order::Integer
    a_q::Matrix{T}
    a_p::Matrix{T}
    α_q::Matrix{T}
    α_p::Matrix{T}

    a_qᵠ::Matrix{T}
    a_pᵠ::Matrix{T}
    α_qᵠ::Matrix{T}
    α_pᵠ::Matrix{T}

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

    function TableauSPARK(order, a_q, a_p, α_q, α_p, a_qᵠ, a_pᵠ, α_qᵠ, α_pᵠ,
                                 b_q, b_p, β_q, β_p,
                                 c_q, c_p, c_λ,
                                 ω_q, ω_p, ω_λ)
        # TODO Make ω_q, ω_p, ω_λ optional arguments.
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(R, Integer)
        @assert isa(order, Integer)
        @assert S > 0
        @assert R > 0
        @assert S==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert S==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        @assert S==size(α_q,1)==size(α_p,1)==length(β_q)==length(β_p)
        @assert R==size(α_q,2)==size(α_p,2)
        @assert R==length(c_λ)
        @assert R==size(a_qᵠ,1)==size(a_pᵠ,1)
        @assert R==size(α_qᵠ,1)==size(α_qᵠ,2)
        @assert R==size(α_pᵠ,1)==size(α_pᵠ,2)
        @assert S==size(a_qᵠ,2)==size(a_pᵠ,2)
        # TODO Add assertions on ω_q, ω_p, ω_λ to be (S-1)x(S) or (R-1)x(R) if set.
        new(order, a_q, a_p, α_q, α_p, a_qᵠ, a_pᵠ, α_qᵠ, α_pᵠ,
                   b_q, b_p, β_q, β_p,
                   c_q, c_p, c_λ,
                   ω_q, ω_p, ω_λ)
    end
end

# TODO Add external constructor for TableauSPARK.

# TODO function readTableauSPARKFromFile(dir::AbstractString, name::AbstractString)


"TableauGLM: Holds the tableau of a general linear method."
immutable TableauGLM{Name, S, R, T} <: Tableau{Name, T}
    order::Integer
    a::Matrix{T}
    b::Matrix{T}
    u::Matrix{T}
    v::Matrix{T}
    c::Vector{T}

    function TableauGLM(order, a, b, u, v, c)
        @assert T <: Real
        @assert isa(Name, Symbol)
        @assert isa(S, Integer)
        @assert isa(R, Integer)
        @assert isa(order, Integer)
        @assert S > 0
        @assert R > 0
        @assert S==length(c)
        @assert S==size(a,1)==size(a,2)==length(u,1)==length(b,2)
        @assert R==size(v,1)==size(v,2)==length(u,2)==length(b,1)
        new(order, a, b, u, v, c)
    end
end

# TODO Add external constructor for TableauGLM.

# TODO function readTableauGLMFromFile(dir::AbstractString, name::AbstractString)


# TODO Add TableauSGALM. # Super General Additive Linear Method
# TODO Add TableauTSRK.  # Two Step Runge Kutta Method
# TODO Add TableauSplitting.

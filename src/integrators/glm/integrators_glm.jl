
"Holds the tableau of a general linear method."
struct TableauGLM{T} <: AbstractTableau{T}
    name::Symbol
    o::Int
    s::Int
    r::Int

    a::Matrix{T}
    b::Matrix{T}
    u::Matrix{T}
    v::Matrix{T}
    c::Vector{T}

    function TableauGLM{T}(name, o, s, r, a, b, u, v, c) where {T}
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

# TODO function writeTableauToFile{Name, T}(dir::AbstractString, tab::TableauGLM{Name, T})


"""
`TableauIPRK`: Tableau of an Implicit Partitioned Runge-Kutta method
```math
\\begin{align*}
P_{n,i} &= \\dfrac{\\partial L}{\\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \\sum \\limits_{j=1}^{s} a_{ij} \\, V_{n,j} , &
q_{n+1} &= q_{n} + h \\sum \\limits_{i=1}^{s} b_{i} \\, V_{n,i} , \\\\
F_{k,i} &= \\dfrac{\\partial L}{\\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h  \\sum \\limits_{i=1}^{s} \\bar{a}_{ij} \\, F_{n,j} , &
p_{n+1} &= p_{n} + h \\sum \\limits_{i=1}^{s} \\bar{b}_{i} \\, F_{n,i} ,
\\end{align*}
```
usually satisfying the symplecticity conditions
```math
\\begin{align*}
b_{i} \\bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\\bar{b}_i &= b_i .
\\end{align*}
```
"""
struct TableauIPRK{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}
    p::CoefficientsRK{T}

    function TableauIPRK{T}(name, o, q, p) where {T}
        @assert q.s==p.s
        new(name, o, q.s, q, p)
    end
end

function TableauIPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}) where {T}
    TableauIPRK{T}(name, order, q, p)
end

# TODO function readTableauIPRKFromFile(dir::AbstractString, name::AbstractString)

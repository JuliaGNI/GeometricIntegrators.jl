
"""
`TableauVPRK`: Tableau of a Variational Partitioned Runge-Kutta method
```math
\\begin{align*}
P_{n,i} &= \\dfrac{\\partial L}{\\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \\sum \\limits_{j=1}^{s} a_{ij} \\, V_{n,j} , &
q_{n+1} &= q_{n} + h \\sum \\limits_{i=1}^{s} b_{i} \\, V_{n,i} , \\\\
F_{k,i} &= \\dfrac{\\partial L}{\\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h  \\sum \\limits_{i=1}^{s} \\bar{a}_{ij} \\, F_{n,j} - d_i \\lambda , &
p_{n+1} &= p_{n} + h \\sum \\limits_{i=1}^{s} \\bar{b}_{i} \\, F_{n,i} , \\\\
&&
0 &= \\sum \\limits_{i=1}^{s} d_i V_i , &&
\\end{align*}
```
satisfying the symplecticity conditions
```math
\\begin{align*}
b_{i} \\bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\\bar{b}_i &= b_i .
\\end{align*}
```
"""
immutable TableauVPRK{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}
    p::CoefficientsRK{T}

    d::Vector{T}

    function TableauVPRK(name, o, q, p, d)
        @assert q.s == p.s == length(d)
        new(name, o, q.s, q, p, d)
    end
end

function TableauVPRK{T}(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d::Vector{T})
    TableauVPRK{T}(name, order, q, p, d)
end

# TODO function readTableauVPRKFromFile(dir::AbstractString, name::AbstractString)

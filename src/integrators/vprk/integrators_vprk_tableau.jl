@doc raw"""
`TableauVPRK`: Tableau of a Variational Partitioned Runge-Kutta method
```math
\begin{aligned}
P_{n,i} &= \dfrac{\partial L}{\partial v} (Q_{n,i}, V_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= \dfrac{\partial L}{\partial q} (Q_{n,i}, V_{n,i}) , &
P_{n,i} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} - d_i \lambda , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} , \\
&&
0 &= \sum \limits_{i=1}^{s} d_i V_i , &&
\end{aligned}
```
satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```
"""
struct TableauVPRK{T} <: AbstractTableau{T}
    @HeaderTableau

    q::Tableau{T}
    p::Tableau{T}

    R∞::Int

    d::Vector{T}

    function TableauVPRK{T}(name, o, q, p, R∞, d) where {T}
        @assert q.s == p.s == length(d)
        new(name, o, q.s, q, p, R∞, d)
    end

    function TableauVPRK{T}(name, o, q, p, R∞) where {T}
        @assert q.s == p.s
        new(name, o, q.s, q, p, R∞)
    end
end

function TableauVPRK(name::Symbol, order::Int, q::Tableau{T}, p::Tableau{T}, R∞::Int, d::Vector{T}) where {T}
    TableauVPRK{T}(name, order, q, p, R∞, d)
end

function TableauVPRK(name::Symbol, order::Int, q::Tableau{T}, p::Tableau{T}, R∞::Int) where {T}
    TableauVPRK{T}(name, order, q, p, R∞)
end

function TableauVPRK(name::Symbol, order::Int, q::Tableau{T}, R∞::Int, d::Vector{T}) where {T}
    TableauVPRK{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞, d)
end

function TableauVPRK(name::Symbol, order::Int, q::Tableau{T}, R∞::Int) where {T}
    TableauVPRK{T}(name, order, q, get_symplectic_conjugate_coefficients(q), R∞)
end

Base.:(==)(tab1::TableauVPRK, tab2::TableauVPRK) = (tab1.o == tab2.o
                                                 && tab1.s == tab2.s
                                                 && tab1.q == tab2.q
                                                 && tab1.p == tab2.p
                                                 && tab1.R∞ == tab2.R∞)

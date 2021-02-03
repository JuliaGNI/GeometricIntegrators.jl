
using LinearAlgebra: istril
using RungeKutta: isexplicit, isimplicit, isdiagnonallyimplicit, isfullyimplicit, istrilstrict

@doc raw"""
Tableau of a Partitioned Runge-Kutta method
```math
\begin{aligned}
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} .
\end{aligned}
```

Parameters:
 * `T`: datatype of coefficient arrays

Fields:
 * `name`: symbolic name of the tableau
 * `o`: order of the method
 * `s`: number of stages
 * `q`: Tableau for `q`
 * `p`: Tableau for `p`

The actual tableaus are stored in `q` and `p`:
 * `a`: coefficients $a_{ij}$ with $ 1 \le i,j \le s$
 * `b`: weights $b_{i}$  with $ 1 \le i \le s$
 * `c`: nodes $c_{i}$  with $ 1 \le i \le s$

Constructors:
```julia
TableauPRK{T}(name, o, s, q, p)
TableauPRK{T}(name, q, p)
TableauPRK(name::Symbol, q::Tableau, p::Tableau)
TableauPRK(name::Symbol, q::Tableau)
```
"""
struct TableauPRK{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::Tableau{T}
    p::Tableau{T}

    function TableauPRK{T}(name, o, s, q, p) where {T}
        @assert s = q.s == p.s
        new(name, o, s, q, p)
    end

    function TableauPRK{T}(name, q, p) where {T}
        new(name, min(q.o, p.o), q.s, q, p)
    end
end

TableauPRK(name::Symbol, q::Tableau{T}, p::Tableau{T}) where {T} = TableauPRK{T}(name, q, p)
TableauPRK(name::Symbol, q::Tableau{T}) where {T} = TableauPRK{T}(name, q, q)

RungeKutta.isexplicit(tab::TableauPRK) = istril(tab.q.a) && istril(tab.p.a) && all([tab.q.a[i,i] == 0 || tab.p.a[i,i] == 0 for i in 1:tab.s]) && (tab.q.c[1] == 0 || tab.p.c[1] == 0)
RungeKutta.isimplicit(tab::TableauPRK) = !isexplicit(tab)
RungeKutta.isdiagnonallyimplicit(tab::TableauPRK) = isimplicit(tab) && tab.s != 1 && istril(tab.q.a) && istril(tab.p.a)
RungeKutta.isfullyimplicit(tab::TableauPRK) = !isexplicit(tab) && !(isdiagnonallyimplicit(tab))


"Holds the tableau of a Runge-Kutta method."
abstract type AbstractTableauRK{T <: Real} <: AbstractTableau{T} end

"Holds the tableau of a partitioned Runge-Kutta method."
abstract type AbstractTableauPRK{T} <: AbstractTableauRK{T} end

Base.:(==)(tab1::AbstractTableauPRK, tab2::AbstractTableauPRK) = (tab1.q == tab2.q && tab1.p == tab2.p && tab1.o == tab2.o && tab1.s == tab2.s)

Base.isequal(tab1::AbstractTableauRK{DT1}, tab2::AbstractTableauRK{DT2}) where {DT1, DT2} = (tab1 == tab2 && DT1 == DT2 && typeof(tab1) == typeof(tab2))

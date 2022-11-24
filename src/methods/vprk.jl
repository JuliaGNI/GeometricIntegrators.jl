
# Variational Partitioned Runge-Kutta Methods

abstract type VPRKMethod <: IODEMethod end


struct VPRK{TT,DT} <: VPRKMethod
    tableau::TT
    d::DT

    function VPRK(tableau::TT, d::DT=nothing) where {TT <: PartitionedTableau, DT <: Union{AbstractVector,Nothing}}
        new{TT,DT}(tableau, d)
    end
end

VPRK(tableau::Tableau, args...; kwargs...) = VPRK(PartitionedTableau(tableau), args...; kwargs...)

tableau(method::VPRK) = method.tableau
nullvector(method::VPRK) = method.d

hasnullvector(method::VPRK{DT,Nothing}) where {DT} = false
hasnullvector(method::VPRK{DT,<:AbstractVector}) where {DT} = true

Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::VPRK; kwargs...) = IntegratorVPRK(problem, tableau(method), nullvector(method); kwargs...)



@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
the symmetric SRK3 coefficients with 3 stages
for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPSRK3(args...; kwargs...) = VPRK(TableauSRK3(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Gauss-Legendre coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKGauss(args...; kwargs...) = VPRK(TableauGauss(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Radau IIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKRadauIIA(args...; kwargs...) = VPRK(TableauRadauIIA(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Radau IIB for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKRadauIIB(args...; kwargs...) = VPRK(TableauRadauIIB(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIA(args...; kwargs...) = VPRK(TableauLobattoIIIA(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIB for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIB(args...; kwargs...) = VPRK(TableauLobattoIIIB(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIC for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIC(args...; kwargs...) = VPRK(TableauLobattoIIIC(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIC̄ for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIC̄(args...; kwargs...) = VPRK(TableauLobattoIIIC̄(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIID for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIID(args...; kwargs...) = VPRK(TableauLobattoIIID(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIE for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIE(args...; kwargs...) = VPRK(TableauLobattoIIIE(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIF for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIF(args...; kwargs...) = VPRK(TableauLobattoIIIF(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIG for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
VPRKLobattoIIIG(args...; kwargs...) = VPRK(TableauLobattoIIIG(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIA for the coefficients $a_{ij}$ and Lobatto IIIB for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIAIIIB(args...; kwargs...) = VPRK(PartitionedTableau(:LobIIIAIIIB, TableauLobattoIIIA(args...), TableauLobattoIIIB(args...)), get_lobatto_nullvector(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIB for the coefficients $a_{ij}$ and Lobatto IIIA for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIBIIIA(args...; kwargs...) = VPRK(PartitionedTableau(:LobIIIBIIIA, TableauLobattoIIIB(args...), TableauLobattoIIIA(args...)), get_lobatto_nullvector(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIA for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIAIIIĀ(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIA(args...)), get_lobatto_nullvector(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIB for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIBIIIB̄(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIB(args...)), get_lobatto_nullvector(args...); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIC for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIICIIIC̄(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIC(args...)); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIC̄ for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIC̄IIIC(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIC̄(args...)); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIID for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIDIIID̄(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIID(args...)); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIE for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIEIIIĒ(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIE(args...)); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIF for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIFIIIF̄(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIF(args...)); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIF̄ for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIF̄IIIF(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIF̄(args...)); kwargs...)


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Lobatto IIIG for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
VPRKLobattoIIIGIIIḠ(args...; kwargs...) = VPRK(SymplecticPartitionedTableau(TableauLobattoIIIG(args...)); kwargs...)

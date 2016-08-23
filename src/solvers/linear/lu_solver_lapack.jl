
import Base.LinAlg: Factorization, checksquare, chkstride1, @assertnonsingular
import Base.LinAlg.BLAS: BlasFloat, BlasInt, liblapack, @blasfunc
import Base.LinAlg.LAPACK: chkargsok, chklapackerror, chktrans


immutable LUSolver{T<:BlasFloat, S<:AbstractMatrix} <: LinearSolver
    factors::S
    pivots::Vector{BlasInt}
    info::Ref{BlasInt}
    LUSolver(factors::AbstractMatrix{T}, ipiv::Vector{BlasInt}, info::Ref{BlasInt}) = new(factors, ipiv, info)
end

function LUSolver(A::AbstractMatrix)
    n = checksquare(A)
    LUSolver{eltype(A), typeof(A)}(A, zeros(BlasInt, n), Ref{BlasInt}())
end

## LAPACK LU factorization and solver for general matrices (GE)
for (getrf, getrs, elty) in
    ((:dgetrf_,:dgetrs_,:Float64),
     (:sgetrf_,:sgetrs_,:Float32),
     (:zgetrf_,:zgetrs_,:Complex128),
     (:cgetrf_,:cgetrs_,:Complex64))
    @eval begin
        # SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, M, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * )
        function getrf!(A::StridedMatrix{$elty}, ipiv::StridedVector{BlasInt}, info::Ref{BlasInt})
            chkstride1(A)
            m = checksquare(A)
            n = m
            lda = m
            # ipiv = similar(A, BlasInt, n)
            # info = Ref{BlasInt}()
            ccall((@blasfunc($getrf), liblapack), Void,
                  (Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
                   &m, &n, A, &lda, ipiv, info)
            chkargsok(info[])
        end

        #     SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          TRANS
        #      INTEGER            INFO, LDA, LDB, N, NRHS
        #     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function getrs!(A::StridedMatrix{$elty}, ipiv::StridedVector{BlasInt}, info::Ref{BlasInt}, b::StridedVector{$elty})
            trans = 'N'
            # chktrans(trans)
            chkstride1(A, b, ipiv)
            n = checksquare(A)
            lda = n
            ldb = length(b)
            if n != ldb
                throw(DimensionMismatch("b has length $ldb, but needs $n"))
            end
            nrhs = 1
            # info = Ref{BlasInt}()
            ccall((@blasfunc($getrs), liblapack), Void,
                  (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}),
                   &trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, info)
            chklapackerror(info[])
        end
    end
end


"""
    getrf!(A, ipiv, info)
Compute the pivoted `LU` factorization of `A`, `A = LU`.
Returns `A`, modified in-place, `ipiv`, the pivoting information, and an `info`
code which indicates success (`info = 0`), a singular value in `U`
(`info = i`, in which case `U[i,i]` is singular), or an error code (`info < 0`).
"""
getrf!(A::StridedMatrix, ipiv::StridedVector{BlasInt}, info::Ref{BlasInt})

"""
    getrs!(A, ipiv, info, b)
Solves the linear equation `A * x = b` for square `A`.
Modifies the vector `b` in place with the solution x. `A` is the `LU`
factorization from `getrf!`, with `ipiv` the pivoting information.
"""
getrs!(A::StridedMatrix, ipiv::StridedVector{BlasInt}, info::Ref{BlasInt}, b::StridedVector)

function lufactorize!{T<:BlasFloat, S<:StridedMatrix}(A::LUSolver{T, S})
    getrf!(A.factors, A.pivots, A.info)
    if A.info[] > 0
        throw(SingularException(A.info[]))
    elseif A.info[] < 0
        throw(ErrorException(A.info[]))
    end
end

function lusolve!{T<:BlasFloat, S<:StridedMatrix}(A::LUSolver{T, S}, b::StridedVector{T})
    getrs!(A.factors, A.pivots, A.info, b)
    if A.info[] < 0
        throw(ErrorException(A.info[]))
    end
end

"""
    AbstractHermite{T,M<:AbstractMatrix{T}} <: Factorization{T}

Supertype for the result of Hermite normal form calculations.
"""
abstract type AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
end

"""
    RowHermite{T<:Integer,M<:AbstractMatrix{T}}

Describes the result of a decomposition of an integer matrix `A` into its row-style Hermite normal
form `H`, which is upper triangular, and a unimodular matrix `U`, such that `A = U/H`.
"""
struct RowHermite{T,M} <: AbstractHermite{T,M}
    H::UpperTriangular{T,M}
    U::M
    info::LinearAlgebra.BlasInt
    function RowHermite(
        H::AbstractMatrix{T1},
        U::M,
        info::Integer = 0
    ) where {T1<:Integer,M<:AbstractMatrix{T2} where T2<:Integer}
        @assert isunimodular(U) "U is not unimodular."
        T = promote_type(T1,T2)
        return new{T,M}(UpperTriangular(H), U, info)
    end
end

"""
    ColumnHermite{T<:Integer,M<:AbstractMatrix{T}}

Describes the result of a decomposition of an integer matrix `A` into its row-style Hermite normal
form `H`, which is lower triangular, and a unimodular matrix `U`, such that `A = U\\H`.
"""
struct ColumnHermite{T,M} <: AbstractHermite{T,M}
    H::LowerTriangular{T,M}
    U::M
    info::LinearAlgebra.BlasInt
    function RowHermite(
        H::AbstractMatrix{T1},
        U::M,
        info::Integer = 0
    ) where {T1<:Integer,M<:AbstractMatrix{T2} where T2<:Integer}
        @assert isunimodular(U) "U is not unimodular."
        T = promote_type(T1,T2)
        return new{T,M}(LowerTriangular(H), U, info)
    end
end

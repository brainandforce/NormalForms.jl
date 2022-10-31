"""
    Smith{T,M<:AbstractMatrix{T}} <: Factorization{T}

Describes the result of the decomposition of an integer matrix `A` into its Smith normal form `S`,
which is diagonal, and unimodular matrices `U` and `V` such that `A = U\\S/V`.
"""
struct Smith{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
    S::Diagonal{T}
    U::M
    V::M
    info::LinearAlgebra.BlasInt
    function Smith(
        S::AbstractMatrix{T1},
        U::AbstractMatrix{T2},
        V::AbstractMatrix{T3},
        info::Integer
    ) where {T1<:Integer,T2<:Integer,T3<:Integer}
        return new{promote_type(T1, T2, T3),typeof(S)}(S, U, V, info)
    end
end

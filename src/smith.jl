"""
    Smith{T,M<:AbstractMatrix{T}} <: Factorization{T}

Describes the result of the decomposition of an integer matrix `A` into its Smith normal form `S`,
which is diagonal, and unimodular matrices `U` and `V` such that `S == U*A*V`, or equivalently,
`A == U\\S/V`.
"""
struct Smith{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
    S::M
    U::M
    V::M
    info::LinearAlgebra.BlasInt
    function Smith(
        S::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        V::AbstractMatrix{<:Integer},
        info::Integer = 0
    )
        @assert isdiag(S) "S is not diagonal."
        @assert isunimodular(U) "U is not unimodular."
        @assert isunimodular(V) "V is not unimodular."
        @assert size(S,1) == size(U,2) string(
            "S and U have incompatible dimensions: ",
            join(size(S), '×'), " for S and ", join(size(U), '×'), " for U"
        )
        @assert size(S,2) == size(V,1) string(
            "S and V have incompatible dimensions: ",
            join(size(S), '×'), " for S and ", join(size(V), '×'), " for V"
        )
        T = Base.promote_eltypeof(S, U, V)
        M = Base.promote_typeof(S, U, V)
        return new{T,M}(S, U, V, info)
    end
end

# Component destructuring
Base.iterate(F::Smith) = (F.S, Val{:U}())
Base.iterate(F::Smith, ::Val{:S}) = iterate(F)
Base.iterate(F::Smith, ::Val{:U}) = (F.U, Val{:V}())
Base.iterate(F::Smith, ::Val{:V}) = (F.V, nothing)
Base.iterate(F::Smith, ::Nothing) = nothing

issuccess(F::Smith) = iszero(F.info)

function Base.summary(io::IO, F::Smith)
    print(io, join(size(F.S), '×'), ' ', typeof(F), ":")
end

function Base.show(io::IO, mime::MIME"text/plain", F::Smith)
    if issuccess(F)
        summary(io, F)
        println(io, "\nHermite normal form:")
        Base.show(io, mime, F.S)
        println(io, "\nLeft unimodular factor:")
        Base.show(io, mime, F.U)
        println(io, "\nRight unimodular factor:")
        Base.show(io, mime, F.V)
    else
        print(io, "Failed factorization of type $(typeof(F))")
    end
end

# TODO: This is not the most efficient way of doing things, nor is this guaranteed to be correct.
"""
    snf!(M::AbstractMatrix{<:Integer}) -> Smith{eltype(M),M}

Calculates the Smith normal form of an integer matrix in-place. This is done naively by applying
`NormalForms.hnf_ma!()` twice, on either side.

!!! warning This may not always work! This is a naive stopgap implementation...
"""
function snf!(M::AbstractMatrix{<:Integer})
    # Just get the unimodular part of the HNF factorization of M
    V = hnf_ma!(M)[2]
    # Now get the unimodular part of the oppposite HNF factorization
    U = hnf_ma!(M')[2]
    # Return as a Smith normal form
    return Smith(M, U, V, 0)
end

"""
    snf(M::AbstractMatrix{<:Integer}) -> Smith{eltype(M),M}

Calculates the Smith normal form of an integer matrix, returning a copy of the original matrix.
This is done naively by applying the Kannan-Bachem algorithm twice, on either side. 

!!! warning This may not always work! This is a naive stopgap implementation...
"""
snf(M::AbstractMatrix{<:Integer}) = snf!(deepcopy(M))

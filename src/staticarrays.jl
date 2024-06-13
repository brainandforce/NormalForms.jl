using .StaticArrays

# Matrix appears slightly faster, but MMatrix allocates less memory
# TODO: perhaps extend this for the sake of StaticArrays?
detb(M::StaticMatrix) = detb!(convert(MMatrix, M))

#= TODO:
    The RowHermite, ColumnHermite, and Smith types assume that the types of the normal form matrix 
    and the types of the U/V factors are identical.

    This is fine for Julia's Matrix type, but it is a serious problem for any StaticMatrix types
    that do not represent square matrices: the U/V factors are going to be square matrices of
    different sizes.

    The stopgap solution is to return `Matrix` types for factorizations of non-square matrices.

    In the future, there are two possible solutions to this issue:
      * Increase the number of type parameters for the factorization types so everything can be
        represented individually.
      * Create specific factorization types for interoperability with StaticArrays.
    I'm inclined to go with the latter (-BF)
=#

# Square matrix methods
function hnfc(M::T, R::RoundingMode = NegativeOffDiagonal) where {D,T<:StaticMatrix{D,D,<:Integer}}
    (H, U, info) = hnf_ma!(MMatrix(M), R)
    return ColumnHermite(T(H), T(U), info)
end

function hnfr(M::T, R::RoundingMode = NegativeOffDiagonal) where {D,T<:StaticMatrix{D,D,<:Integer}}
    (H, U, info) = hnf_ma!(MMatrix(M'), R)
    return ColumnHermite(T(H'), T(U'), info)
end

function snf(M::T) where {D,T<:StaticMatrix{D,D,<:Integer}}
    (S, U, V, info) = snf_ma!(MMatrix(M))
    return Smith(T(S), T(U), T(V), info)
end

# Non-square matrix methods
function hnfc(M::StaticMatrix{<:Any,<:Any,<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    return hnfc!(collect(M), R)
end

function hnfr(M::StaticMatrix{<:Any,<:Any,<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    return hnfr!(collect(M), R)
end

snf(M::StaticMatrix{<:Any,<:Any,<:Integer}) = snf!(collect(M))

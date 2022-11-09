"""
    AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}

Supertype for the result of Hermite normal form calculations.
"""
abstract type AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
end

# Component destructuring
Base.iterate(F::AbstractHermite) = (F.H, Val{:U}())
Base.iterate(F::AbstractHermite, ::Val{:H}) = iterate(F)
Base.iterate(F::AbstractHermite, ::Val{:U}) = (F.U, nothing)
Base.iterate(F::AbstractHermite, ::Nothing) = nothing

LinearAlgebra.issuccess(H::AbstractHermite) = iszero(H.info)
LinearAlgebra.diag(H::AbstractHermite) = diag(H.H)

"""
    RowHermite{T<:Integer,M<:AbstractMatrix{T}}

Describes the result of a decomposition of an integer matrix `A` into its row-style Hermite normal
form `H`, which is upper triangular, and a unimodular matrix `U`, such that `H == U*A`, or
equivalently, `A = U\\H`.
"""
struct RowHermite{T,M} <: AbstractHermite{T,M}
    H::M
    U::M
    info::LinearAlgebra.BlasInt
    function RowHermite(
        H::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        info::Integer = 0
    )
        @assert istriu(H) "H is not upper triangular."
        @assert isunimodular(U) "U is not unimodular."
        @assert size(H,1) == size(U,1) string(
            "H and U have incompatible dimensions: ",
            join(size(H), '×'), " for H and ", join(size(U), '×'), "for U"
        )
        M = promote_typeof(H,U)
        return new{eltype(M),M}(promote(H, U)..., info)
    end
end

"""
    ColumnHermite{T<:Integer,M<:AbstractMatrix{T}}

Describes the result of a decomposition of an integer matrix `A` into its row-style Hermite normal
form `H`, which is lower triangular, and a unimodular matrix `U`, such that `H == A*U`, or
equivalently, `A == H/U`.
"""
struct ColumnHermite{T,M} <: AbstractHermite{T,M}
    H::M
    U::M
    info::LinearAlgebra.BlasInt
    function ColumnHermite(
        H::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        info::Integer = 0
    )
        @assert istril(H) "H is not lower triangular."
        @assert isunimodular(U) "U is not unimodular."
        @assert size(H,2) == size(U,1) string(
            "H and U have incompatible dimensions: ",
            join(size(H), '×'), " for H and ", join(size(U), '×'), "for U"
        )
        M = promote_typeof(H, U)
        return new{eltype(M),M}(promote(H, U)..., info)
    end
end

function Base.summary(io::IO, H::AbstractHermite)
    print(io, join(size(H.H), '×'), ' ', typeof(H), ":")
end

function Base.show(io::IO, mime::MIME"text/plain", F::AbstractHermite)
    if issuccess(F)
        summary(io, F)
        println(io, "\nHermite normal form:")
        Base.show(io, mime, F.H)
        println(io, "\nUnimodular factor:")
        Base.show(io, mime, F.U)
    else
        print(io, "Failed factorization of type $(typeof(F))")
    end
end

"""
    PositiveOffDiagonal

Alias for `RoundDown`. In `NormalForms.reduce_cols_off_diagonal!()`, this makes all off-diagonal
elements the smallest positive value possible.
"""
const PositiveOffDiagonal = RoundDown

"""
    NegativeOffDiagonal

Alias for `RoundUp`. In `NormalForms.reduce_cols_off_diagonal!()`, this makes all off-diagonal
elements the smallest positive value possible.
"""
const NegativeOffDiagonal = RoundUp

"""
    NormalForms._hnf_ma!(
        A::AbstractMatrix{<:Integer},
        R::RoundingMode = NegativeOffDiagonal
    ) -> Tuple{typeof(A), typeof(A), LinearAlgebra.BlasInt}

A modified version of `HermiteNormalForm._hnf_like!()` that calculates the column-style Hermite
normal form of a matrix in-place, returning `A`, the unimodular matrix `U` such that `H == A*U`,
and the factorization info (zero if success - matching LinearAlgebra factorizations).

The original implementation of `HermiteNormalForm._hnf_like!()` was written by Yingbo Ma, and may
be found here: https://github.com/YingboMa/HermiteNormalForm.jl/
"""
function hnf_ma!(
    A::AbstractMatrix{<:Integer},
    R::RoundingMode = NegativeOffDiagonal
)
    # Create the unimodular matrix as an identity matrix
    U = eye(A,2)
    # Convert to a transpose 
    if A isa Adjoint
        U = U'
    elseif A isa Transpose
        U = transpose(U)
    end
    # Loop through each row of A
    @inbounds for k in axes(A,1)
        pivot = find_pivot(A,k)
        # It's not entirely clear to me what ki does.
        if iszero(pivot)
            # Set ki to the first nonzero column member
            ki = findfirst(!iszero, @view A[k+1:end, k])
            # If all are zero, break from the row loop (matrix is truly rank deficient)
            # Otherwise, add the found value to k to compensate for the view's indices
            isnothing(ki) ? continue : (ki += k) # what does it mean when ki > k...
        else
            ki = k
            swapcols!(A, pivot, k)
            swapcols!(U, pivot, k)
        end
        zero_row!(A, U, k, ki)
        @debug string(
            "Matrix state after zeroing along row or column k = $k:",
            "\nA = ", repr("text/plain", A),
            "\nU = ", repr("text/plain", A),
        )
        @assert isunimodular(U)
        # Skip the extra steps if k is not in the largest leading minor
        k <= size(A,2) && reduce_cols_off_diagonal!(A, U, k, ki, R)
        @debug string(
            "Matrix state after reducing off-diagonal at k = $k:",
            "\nA = ", repr("text/plain", A),
            "\nU = ", repr("text/plain", A),
        )
        @assert isunimodular(U)
    end
    return (A, U, 0)
end

"""
    hnfc!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> ColumnHermite{eltype(M),typeof(M)}

Calculates the column Hermite normal form of the integer matrix `M` in-place, returning a
`ColumnHermite` describing the elements of the factorization.
"""
function hnfc!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    return ColumnHermite(hnf_ma!(M, R)...)
end

"""
    hnfc(M::AbstractMatrix{T<:Integer}) -> ColumnHermite{T,typeof(M)}
        -> ColumnHermite{eltype(M),typeof(M)}

Calculates the column Hermite normal form of the integer matrix `M` in-place, returning a
`ColumnHermite` describing the elements of the factorization. Unlike `hnfc!()`, this function
creates a copy of the input matrix rather than modifying it in-place.
"""
hnfc(M::AbstractMatrix{<:Integer}) = hnfc!(deepcopy(M))

"""
    hnfr!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> RowHermite{eltype(M),typeof(M)}

Calculates the row Hermite normal form of the integer matrix `M` in-place, returning a `RowHermite`
describing the elements of the factorization.
"""
function hnfr!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    # Calculate the column HNF of the transpose
    (H, U, info) = hnf_ma!(M', R)
    # Transpose back to get M again
    return RowHermite(H', U', info)
end

"""
    hnfr(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> RowHermite{eltype(M),typeof(M)}

Calculates the row Hermite normal form of the integer matrix `M` in-place, returning a `RowHermite`
describing the elements of the factorization. Unlike `hnfr!()`, this function creates a copy of the
input matrix rather than modifying it in-place.
"""
hnfr(M::AbstractMatrix{<:Integer}) = hnfr!(deepcopy(M))

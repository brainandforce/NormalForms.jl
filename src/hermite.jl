"""
    AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}

Supertype for the result of Hermite normal form calculations.
"""
abstract type AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
end

# Component destructuring
Base.iterate(F::AbstractHermite) = (F.H, Val{:U}())
Base.iterate(F::AbstractHermite, ::Val{:H}) = iterate(F)
Base.iterate(F::AbstractHermite, ::Val{:U}) = (F.H, nothing)
Base.iterate(F::AbstractHermite, ::Nothing) = nothing

issuccess(H::AbstractHermite) = iszero(H.info)

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
        M = promote_type(typeof(H), typeof(U))
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
        M = Base.promote_typeof(H, U)
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

const PositiveOffDiagonal = RoundDown
const NegativeOffDiagonal = RoundUp

#=
"""
    NormalForms.kb_bezout(x::Integer, y::Integer)

The modified Euclidean algorithm used in the Kannan-Bachem algorithm for generating the elementary
matrix operations used in the main loop.

# References

(1) Kannan, R.; Bachem, A. Polynomial Algorithms for Computing the Smith and Hermite Normal Forms
of an Integer Matrix. SIAM J. Comput. 1979, 8 (4), 499–507. https://doi.org/10.1137/0208040.
"""
function kb_bezout(x::Integer, y::Integer)
    # Order a and b as needed
    (a,b) = (x,y) .|> (minimum, maximum)
    (r,p,q) = gcdx(a,b)
    # Modify p and q if needed
    if abs(q) > abs(a)
        p += fld(q,a) * b
        q -= fld(q,a) * a
    end
    return (r,p,q)
end
=#

"""
    HermiteStyle{S}

Determines whether the output of `NormalForms.hnf_kb!()` produces a `RowHermite` or `ColumnHermite`
output. The possible choices are `RowStyle` (`HermiteStyle{:row}()`) or `ColumnStyle`
(`HermiteStyle{:col}()`).
"""
struct HermiteStyle{S}
end

"""
    NormalForms._hnf_ma!(
        A::AbstractMatrix{<:Integer},
        R::RoundingMode = NegativeOffDiagonal,
        ::Val{diagonalize} = Val{true}(),
    ) -> Tuple{typeof(A), typeof(A), LinearAlgebra.BlasInt}

A modified version of `HermiteNormalForm._hnf_like!()` that calculates the column-style Hermite
normal form of a matrix in-place, returning `A`, the unimodular matrix `U` such that `H == AU`, and
the factorization info (zero if success - matching LinearAlgebra factorizations).

The original implementation may be found here: https://github.com/YingboMa/HermiteNormalForm.jl/
"""
function hnf_ma!(
    A::AbstractMatrix{<:Integer},
    R::RoundingMode = NegativeOffDiagonal,
    ::Val{diagonalize} = Val{false}(),
) where {diagonalize}
    # Create the unimodular matrix as an identity matrix
    U = diagm(ones(eltype(A), size(A,2)))
    # Convert to a transpose 
    if A isa Union{Adjoint,Transpose}
        U = U'
    end
    # Loop through each row of A
    @inbounds for k in axes(A,1)
        # Set k as a default pivot
        pivot = k
        # If k has a diagonal and it's zero, change the pivot
        if k <= minimum(size(A)) && iszero(A[k,k])
            # Loop through each element of row k
            for j in axes(A,2)[k+1:end]
                # Set the pivot to a nonzero element
                if A[k,j] != zero(eltype(A))
                    pivot = j
                    # We probably should break here to get the *first* nonzero element...
                end
            end
            # If no pivot is found, set the pivot to 0 to indicate rank deficiency
            pivot == k && (pivot = zero(pivot))
        end
        # It's not entirely clear to me what ki does.
        if iszero(pivot)
            # Set ki to the first nonzero column member
            ki = findfirst(!iszero, @view A[k+1:end, k])
            # If all are zero, break from the row loop (matrix is truly rank deficient)
            # Otherwise, add the found value to k to compensate for the view's indices
            isnothing(ki) ? continue : (ki += k) # what does it mean when ki > k...
        else
            ki = k
            # Permute rows according to the pivot if needed
            if pivot != k
                A[:,[pivot,k]] = A[:,[k,pivot]]
                U[:,[pivot,k]] = U[:,[k,pivot]]
            end
        end
        # Zero the off-diagonal elements: A[k,1:k-1]
        for j in axes(A,2)[k+1:end]
            # Why do we use ki here?
            Akk, Akj = A[ki, k], A[ki, j]
            (d, p, q) = gcdx(Akk, Akj)
            Akkd, Akjd = div(Akk, d), div(Akj, d)
            # Mutating A
            Ak, Aj = A[:,k], A[:,j]
            A[:,k] =  Ak * p + Aj * q
            A[:,j] = -Ak * Akjd + Aj * Akkd
            # Mutating U
            Uk, Uj = U[:,k], U[:,j]
            U[:,k] =  Uk * p + Uj * q
            U[:,j] = -Uk * Akjd + Uj * Akkd
        end
        # Skip the extra steps if k is not in the largest leading minor
        k <= size(A,2) || continue
        # Does this calculate a different transform?
        if diagonalize
            # Zero out A[k, 1:k-1] === A21[1, 1:k-1] by doing
            # A[:, j] = A[:, [k j]] * [-A[k, j], A[k, k]]
            #
            # Note that we then have:
            #   A[k, j] = A[k, [k j]] * [-A[k, j], A[k, k]]
            # = A[k, k] * -A[k, j] + A[k, j] * A[k, k] = 0
            for j in 1:k-1
                Akk, Akj = A[ki,k], A[ki, j]
                d = gcd(A[ki,k], A[ki, j])
                Akkd, Akjd = div(A[ki,k], d), div(A[ki,j], d)

                Ak, Aj = A[:,k], A[:,j]
                A[:,j] = -Ak * Akjd + Aj * Akkd

                Uk, Uj = U[:,k], U[:,j]
                U[:,j] = -Uk * Akjd + Uj * Akkd
            end
        else
            # Make the diagonal element positive
            if A[ki, k] < zero(eltype(A))
                @. A[:,k] = -A[:,k]
                @. U[:,k] = -U[:,k]
            end
            # Minimize the off-diagonal elements
            for j = 1:k-1
                mul = div(A[ki,j], A[ki,k], R)
                @. A[:,j] -= mul * A[:,k]
                @. U[:,j] -= mul * U[:,k]
            end
        end
    end
    return (A, U, 0)
end

"""
    NormalForms.hnfc!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
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

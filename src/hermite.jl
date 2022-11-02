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
        M = promote_type(typeof(H), typeof(U))
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

# TODO: Verify this works for sparse matrices...
"""
    NormalForms.hnf_kb!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> Tuple{typeof(M),typeof(M),LinearAlgebra.BlasInt}

Generic algorithm for calculating the column Hermite normal form given a matrix `M` and the
optional specification of a rounding mode, either `NegativeOffDiagonal` (alias for RoundUp, 
default) or `PositiveOffDiagonal` (alias for `RoundDown`).

# References

(1) Kannan, R.; Bachem, A. Polynomial Algorithms for Computing the Smith and Hermite Normal Forms
of an Integer Matrix. SIAM J. Comput. 1979, 8 (4), 499–507. https://doi.org/10.1137/0208040.
"""
function hnf_kb!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    # Generate an identity matrix to start
    # This matrix will remain unimodular
    U = diagm(ones(eltype(M), size(M,2)))
    # Off-diagonal reduction of M, with corresponding changes also made to U
    function rod!(k::Integer)
        # Make diagonal elements positive
        if M[k,k] < 0
            M[:,k] *= -1
            U[:,k] *= -1
        end
        # Modify all rows before k
        # If fld (RoundDown) is used, all off-diagonal elements will be positive
        # If cld (RoundUp) is used, they'll be negative
        for z in 1:(k-1)
            M[:,z] -= div(M[k,z], M[k,k], R) * M[:,k]
            U[:,z] -= div(M[k,z], M[k,k], R) * U[:,k]
        end
    end
    # Permutation vector to track pivots
    P = collect(1:size(M,2))
    # Permute M so that all of its principal minors are nonsingular
    # First, check that the first entry of M is nonzero
    if iszero(M[1,1])
        # Permute so that M[1,1] is not zero;
        # Propagate changes to U and P; M was originally M * U
        let i = findfirst(!iszero, M[1,:])
            M[:,[1,i]] = M[:,[i,1]]
            U[:,[1,i]] = U[:,[i,1]]
            P[[1,i]] = P[[i,1]]
        end
    end
    # Now check the determinants of the square minors
    for i in 2:minimum(size(M))
        # Check that the i×i principal minor is not singular
        if iszero(detb(M[1:i,1:i]))
            # @debug "Singular minor: " * repr("text/plain", M[1:i,1:i])
            # If it is singular, loop through each 
            for j in eachindex(eachcol(M))[(i+1):end]
                # Find a column that makes the determinant nonzero
                if !iszero(detb(M[1:i, vcat(1:(i-1), j)]))
                    # Swap columns i and j
                    M[:,[i,j]] = M[:,[j,i]]
                    U[:,[i,j]] = U[:,[j,i]]
                    P[[i,j]] = P[[j,i]]
                else
                    # if one can't be found, the matrix is singular
                    throw(SingularException)
                end
            end
        end
    end
    # Loop through the rows of M
    for i in eachindex(eachcol(M))[2:end]
        # Loop through columns up to (but excluding) i
        for j in 1:min(i-1, size(M,1))
            @debug "Setting M[$j,$i] to zero..."
            (a,b) = ((M[j,j], M[j,i])) .|> (minimum, maximum)
            # Calculate gcd and Bezout coefficients
            (r,p,q) = gcdx(a,b)
            # Modify p and q if needed
            if abs(q) > abs(a)
                p += fld(q,a) * b
                q -= fld(q,a) * a
            end
            # Make M[j,i] zero and update U to match the transform
            # div() will have no remainder since it divides by the gcd
            # So the rounding mode is irrelevant
            D = [p -div(M[j,i], r); q div(M[j,j], r)]
            M[:,[j,i]] = M[:,[j,i]] * D
            U[:,[j,i]] = U[:,[j,i]] * D
            #=
            @debug string(
                "At i = $i, j = $j:",
                "\nM is now ", repr("text/plain", M),
                "\nU is now ", repr("text/plain", U),
                "\nM[$j,$i] is ", "not "^iszero(M[i,j]), "zero"
            )
            =#
            j > 1 && rod!(j)
        end
        i <= minimum(size(M)) && rod!(i)
        # @debug "Matrix state: " * repr("text/plain", M)
    end
    return (M, U, 0)
end

"""
    NormalForms.hnfc!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> ColumnHermite{eltype(M),typeof(M)}

Calculates the column Hermite normal form of the integer matrix `M` in-place using an improved
Kannan-Bachem algorithm, returning a `ColumnHermite` describing the elements of the factorization.

# References

(1) Kannan, R.; Bachem, A. Polynomial Algorithms for Computing the Smith and Hermite Normal Forms
of an Integer Matrix. SIAM J. Comput. 1979, 8 (4), 499–507. https://doi.org/10.1137/0208040.
"""
function hnfc!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    return ColumnHermite(hnf_kb!(M, R)...)
end

"""
    hnfc(M::AbstractMatrix{T<:Integer}) -> ColumnHermite{T,typeof(M)}
        -> ColumnHermite{eltype(M),typeof(M)}

Calculates the column Hermite normal form of the integer matrix `M` using an improved Kannan-Bachem 
algorithm, returning a `ColumnHermite` describing the elements of the factorization. Unlike
`hnfc!()`, this creates a copy of the input matrix rather than modifying it in-place.

# References

(1) Kannan, R.; Bachem, A. Polynomial Algorithms for Computing the Smith and Hermite Normal Forms
of an Integer Matrix. SIAM J. Comput. 1979, 8 (4), 499–507. https://doi.org/10.1137/0208040.
"""
hnfc(M::AbstractMatrix{<:Integer}) = hnfc!(deepcopy(M))

"""
    hnfr!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> RowHermite{eltype(M),typeof(M)}

Calculates the row Hermite normal form of the integer matrix `M` in-place using an improved
Kannan-Bachem algorithm, returning a `RowHermite` describing the elements of the factorization.

# References

(1) Kannan, R.; Bachem, A. Polynomial Algorithms for Computing the Smith and Hermite Normal Forms
of an Integer Matrix. SIAM J. Comput. 1979, 8 (4), 499–507. https://doi.org/10.1137/0208040.
"""
function hnfr!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
    # Calculate the column Hermite normal form of the transpose
    (H, U, info) = hnf_kb!(M', R)
    # Transpose H back (equal to mutated M)
    return RowHermite(H', U, info)
end

"""
    hnfr(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)
        -> RowHermite{eltype(M),typeof(M)}

Calculates the row Hermite normal form of the integer matrix `M` in-place using an improved
Kannan-Bachem algorithm, returning a `RowHermite` describing the elements of the factorization. 
Unlike `hnfr!()`, this creates a copy of the input matrix rather than modifying it in-place.

# References

(1) Kannan, R.; Bachem, A. Polynomial Algorithms for Computing the Smith and Hermite Normal Forms
of an Integer Matrix. SIAM J. Comput. 1979, 8 (4), 499–507. https://doi.org/10.1137/0208040.
"""
hnfr(M::AbstractMatrix{<:Integer}) = hnfr!(deepcopy(M))

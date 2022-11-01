"""
    AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}

Supertype for the result of Hermite normal form calculations.
"""
abstract type AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
end

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
        H::AbstractMatrix{T1},
        U::AbstractMatrix{T2},
        info::Integer = 0
    ) where {T1<:Integer,T2<:Integer}
        @assert istriu(H) "H is not upper triangular."
        @assert isunimodular(U) "U is not unimodular."
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
        H::AbstractMatrix{T1},
        U::AbstractMatrix{T2},
        info::Integer = 0
    ) where {T1<:Integer,T2<:Integer}
        @assert istril(H) "H is not lower triangular."
        @assert isunimodular(U) "U is not unimodular."
        M = promote_type(typeof(H), typeof(U))
        return new{eltype(M),M}(promote(H, U)..., info)
    end
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
    HermiteStyle{Side,Sign}

A type used for controlling the mode by which a Hermite normal form calculation is calculated.
`Side` determines whether the result is upper triangular (`:Upper`) or lower triangular (`:Lower`).
`Sign` determines whether off-diagonal elements are positive (`:Positive`) or negative
(`:Negative`).

The available styles are:
  * UpperPositive (upper triangular, positive off-diagonal elements)
  * UpperNegative (upper triangular, negative off-diagonal elements) (default)
  * LowerPositive (lower triangular, positive off-diagonal elements)
  * LowerNegative (lower triangular, negative off-diagonal elements)

"""
struct HermiteStyle{Side,Sign}
end

const UpperPositive = HermiteStyle{:Upper,:Positive}()
const UpperNegative = HermiteStyle{:Upper,:Negative}()
const LowerPositive = HermiteStyle{:Lower,:Positive}()
const LowerNegative = HermiteStyle{:Lower,:Negative}()

# TODO: Does the algorithm also work for sparse matrices?
"""
    hnfc!(M::DenseMatrix{T}) -> ColumnHermite{T,typeof(M)}

Calculates the column Hermite normal form of the integer matrix `M` in-place using an improved
Kannan-Bachem algorithm, returning a `ColumnHermite` describing the elements of the factorization.
"""
function hnfc!(M::DenseMatrix{T}) where T<:Integer
    # Generate an identity matrix to start
    # This matrix will remain unimodular
    U = diagm(ones(T, size(M,2)))
    # Off-diagonal reduction of M, with corresponding changes also made to U
    function rod!(k)
        # Make diagonal elements positive
        if M[k,k] < 0
            M[:,k] *= -1
            U[:,k] *= -1
        end
        # Modify all rows before k
        # If fld is used, all off-diagonal elements will be positive
        # If cld is used, they'll be negative
        for z in 1:(k-1)
            M[:,z] -= cld(M[k,z], M[k,k]) * M[:,k]
            U[:,z] -= cld(M[k,z], M[k,k]) * U[:,k]
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
    # Now check the principal minors
    for i in eachindex(eachrow(M))[2:end]
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
    #=
    @debug string(
        "Permuted M to: ", repr("text/plain", M), "\n",
        "Permuted U to: ", repr("text/plain", U), "\n",
        "Permutation vector: ", repr(P)
    )
    =#
    # Save the permutation matrix for later...
    # pmat = deepcopy(U)
    # @info "pmat = $pmat"
    # Loop through the rows of M
    for i in eachindex(eachcol(M))[2:end]
        # Loop through columns up to (but excluding) i
        for j in 1:(i-1)
            @debug "Setting M[$j,$i] to zero..."
            (a,b) = ((M[j,j], M[j,i])) .|> (minimum, maximum)
            # Calculate gcd and Bezout coefficients
            (r,p,q) = gcdx(a,b)
            # Modify p and q if needed
            if abs(q) > abs(a)
                p += fld(q,a) * b
                q -= fld(q,a) * a
            end
            # Make M[j,i] zero
            D = [p -div(M[j,i], r); q div(M[j,j], r)]
            # @debug "D: " * repr("text/plain", D)
            M[:,[j,i]] = M[:,[j,i]] * D
            # Update U to match the transform
            U[:,[j,i]] = U[:,[j,i]] * D
            #=
            @debug string(
                "At i = $i, j = $j:",
                "\nM is now ", repr("text/plain", M),
                "\nU is now ", repr("text/plain", U),
                "\nM[$j,$i] is ", "not "^iszero(M[i,j]), "zero"
            )
            =#
            # Reduce off diagonal
            # @debug "Reducing off diagonal at j = $j"
            j > 1 && rod!(j)
        end
        # Reduce off diagonal once more
        # @debug "Reducing off diagonal at i = $i"
        i < size(M,2) && rod!(i)
        # @debug "Matrix state: " * repr("text/plain", M)
    end
    # Check that the result actually makes sense...
    @assert istril(M) "Result was not lower triangular! Got " * repr("text/plain", M)
    # Use the permutation data to put U in a more convenient form
    return ColumnHermite(M, U, 0)
end

hnfc(M::DenseMatrix{<:Integer}) = hnfc!(deepcopy(M))

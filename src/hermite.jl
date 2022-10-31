"""
    AbstractHermite{T,M<:AbstractMatrix{T}} <: Factorization{T}

Supertype for the result of Hermite normal form calculations.
"""
abstract type AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
end

issuccess(H::AbstractHermite) = H.info == 0

"""
    RowHermite{T<:Integer,M<:AbstractMatrix{T}}

Describes the result of a decomposition of an integer matrix `A` into its row-style Hermite normal
form `H`, which is upper triangular, and a unimodular matrix `U`, such that `A == U/H`.
"""
struct RowHermite{T,M} <: AbstractHermite{T,M}
    H::UpperTriangular{T,M}
    U::M
    info::LinearAlgebra.BlasInt
    function RowHermite(
        H::AbstractMatrix{T0},
        U::M,
        info::Integer = 0
    ) where {T0<:Integer,M<:AbstractMatrix}
        @assert isunimodular(U) "U is not unimodular."
        T = promote_type(T0,eltype(M))
        return new{T,M}(UpperTriangular(H), U, info)
    end
end

"""
    ColumnHermite{T<:Integer,M<:AbstractMatrix{T}}

Describes the result of a decomposition of an integer matrix `A` into its row-style Hermite normal
form `H`, which is lower triangular, and a unimodular matrix `U`, such that `A == U\\H`.
"""
struct ColumnHermite{T,M} <: AbstractHermite{T,M}
    H::LowerTriangular{T,M}
    U::M
    info::LinearAlgebra.BlasInt
    function ColumnHermite(
        H::AbstractMatrix{T0},
        U::M,
        info::Integer = 0
    ) where {T0<:Integer,M<:AbstractMatrix}
        @assert isunimodular(U) "U is not unimodular."
        T = promote_type(T0,eltype(M))
        return new{T,M}(LowerTriangular(H), U, info)
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
    NormalForms.hnfc!(M::DenseMatrix{T}) -> ColumnHermite{T,typeof(M)}

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
        for z in 1:(k-1)
            M[:,z] -= cld(M[k,z], M[k,k]) * M[:,k]
            U[:,z] -= cld(M[k,z], M[k,k]) * U[:,k]
        end
    end
    # Permutation vector to track pivots
    P = collect(1:size(M,2))
    # View into M with permutation vector
    # V = view(M, 1:size(M,1), P)
    # Permute M so that all of its principal minors are nonsingular
    for i in eachindex(eachrow(M))[2:end]
        # Check that the i×i principal minor is not singular
        if iszero(detb(M[1:i,1:i]))
            @debug "Singular minor: " * repr("text/plain", M[1:i,1:i])
            # If it is singular, loop through each 
            for j in eachindex(eachcol(M))[(i+1):end]
                # Find a column that makes the determinant nonzero
                if !iszero(detb(M[1:i, vcat(1:(i-1), j)]))
                    # Swap columns i and j
                    M[:,[i,j]] = M[:,[j,i]]
                    U[:,[i,j]] = U[:,[j,i]]
                    P[[i,j]] = P[[j,i]]
                    @debug string(
                        "Permuted M to: ", repr("text/plain", M), "\n",
                        "Permutation vector: ", repr(P)
                    )
                else
                    # if one can't be found, the matrix is singular
                    throw(SingularException)
                end
            end
        end
    end
    # Loop through the rows of M
    for i in eachindex(eachrow(M))[2:end]
        # Loop through columns up to (but excluding) i
        for j in 1:(i-1)
            @debug "Matrix state: " * repr("text/plain", M)
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
            @debug "D: " * repr("text/plain", D)
            M[:,[j,i]] = M[:,[j,i]] * D
            # Update U to match the transform
            U[:,[j,i]] = U[:,[j,i]] * D
            @debug string(
                "At i = $i, j = $j:",
                "\nM is now ", repr("text/plain", M),
                "\nU is now ", repr("text/plain", U),
                "\nM[$j,$i] is ", "not "^iszero(M[i,j]), "zero"
            )
            # Reduce off diagonal
            j > 1 && rod!(j)
        end
        # Reduce off diagonal once more
        rod!(i)
        @debug "Matrix state: " * repr("text/plain", M)
    end
    @assert istril(M) "Result was not lower triangular! " * repr("text/plain", M)
    # Use the permutation data to put U in a more convenient form
    @info "Permutation vector: $P"
    return ColumnHermite(M, U, 0)
end

hnfc(M::DenseMatrix{<:Integer}) = hnfc!(deepcopy(M))

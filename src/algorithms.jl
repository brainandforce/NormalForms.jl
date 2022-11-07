"""
    NormalForms.eye(T::Type, sz)

Generates an identity matrix of size `sz`.
"""
eye(T::Type, sz) = Matrix{T}(LinearAlgebra.I(sz))

"""
    NormalForms.eye(M::AbstractMatrix{T}, dim)

Generates an identity matrix of the same size as dimension `dim` of `M`.
"""
eye(M::AbstractMatrix{T}, dim) where T = eye(T, size(M)[dim])

"""
    isunimodular(M::AbstractMatrix)

Returns `true` if a matrix is unimodular. A unimodular matrix is a square integer matrix with A
determinant equal to 1 or -1.
"""
isunimodular(M::AbstractMatrix{<:Integer}) = size(M,1) == size(M,2) && isone(abs(detb(M)))

function isunimodular(M::AbstractMatrix)
    try
        isunimodular(Int.(M))
    catch e
        e isa InexactError ? false : throw(e)
    end
end

# TODO: a unimodular inverse?

"""
    NormalForms.gcd_kb(x::Integer, y::Integer) -> NTuple{3,Integer}

Calculates the gcd and Bézout coefficients with the Kannan-Bachem modification.
"""
function gcd_kb(x::Integer, y::Integer)
    (r,p,q) = gcdx(x,y)
    # Check if 
    if q >= x && q >= y
        (a,b) = ((x,y),) .|> (minimum, maximum)
        p = p + fld(q,a) * b
        q = q - fld(q,a) * a
    end
    return (r,p,q)
end

"""
    NormalForms.find_pivot(A::AbstractMatrix{<:Integer}, k::Integer)

Finds a suitable pivot element for the Smith or Hermite normal form algorithms. If none can be
found, it returns 0.
"""
@inline function find_pivot(A::AbstractMatrix{<:Integer}, k::Integer)
    # Default to k as the pivot
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
    return pivot
end

"""
    NormalForms.zero_row!(
        A::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        k::Integer
        ki::Integer = k
    )

Zeroes the elements along the row `A[k,k+1:end]` (all elements right of `A[k,k]`). Changes to `A`
are tracked in the matrix `U`, which will only undergo unimodular transforms.
"""
@inline function zero_row!(
    A::AbstractMatrix{<:Integer},
    U::AbstractMatrix{<:Integer},
    k::Integer,
    ki::Integer = k
)
    # This can safely be skipped if k is the last element
    k < last(axes(A,2)) && for j in axes(A,2)[k+1:end]
        (d, p, q) = gcd_kb(A[ki, k], A[ki, j])
        Akkd, Akjd = div(A[ki, k], d), div(A[ki, j], d)
        # Mutating A
        Ak, Aj = A[:,k], A[:,j]
        A[:,k] =  Ak * p + Aj * q
        A[:,j] = -Ak * Akjd + Aj * Akkd
        # Mutating U
        Uk, Uj = U[:,k], U[:,j]
        U[:,k] =  Uk * p + Uj * q
        U[:,j] = -Uk * Akjd + Uj * Akkd
        # @assert isunimodular(U)
    end
end

"""
    NormalForms.zero_col!(
        A::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        k::Integer
        ki::Integer = k
    )

Zeroes the elements along the column `A[k+1:end,k]` (all elements below `A[k,k]`). Changes to `A`
are tracked in the matrix `U`, which will only undergo unimodular transforms.
"""
@inline function zero_col!(
    A::AbstractMatrix{<:Integer},
    U::AbstractMatrix{<:Integer},
    k::Integer,
    ki::Integer = k
)
    # This can safely be skipped if k is the last element
    k < last(axes(A,1)) && for j in axes(A,1)[k+1:end]
        (d, p, q) = gcd_kb(A[k, ki], A[j, ki])
        Akkd, Ajkd = div(A[k, ki], d), div(A[j, ki], d)
        # Mutating A to zero the upper off-diagonal elements
        Ak, Aj = A[k,:], A[j,:]
        A[k,:] =  Ak * p + Aj * q
        A[j,:] = -Ak * Ajkd + Aj * Akkd
        # Mutating U as the matching unimodular matrix
        Uk, Uj = U[k,:], U[j,:]
        U[k,:] =  Uk * p + Uj * q
        U[j,:] = -Uk * Ajkd + Uj * Akkd
        # @assert isunimodular(U)
    end
end

"""
    NormalForms.reduce_off_diagonal!(
        A::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        R::RoundingMode
    )

Performs an off-diagonal reduction of elements in `A` in-place by columns. This also ensures that
diagonal elements are positive by multiplying the columns by -1 if needed. Changes to `A` are
tracked in the matrix `U`, which will only undergo unimodular transforms.

If the `RoundingMode` is set to `RoundUp` or `NegativeOffDiagonal`, the off-diagonal elements will
be negative. If set to `RoundDown` or `PositiveOffDiagonal`, they will be positive.
"""
@inline function reduce_cols_off_diagonal!(
    A::AbstractMatrix{<:Integer},
    U::AbstractMatrix{<:Integer},
    k::Integer,
    ki::Integer = k,
    R::RoundingMode = NegativeOffDiagonal
)
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

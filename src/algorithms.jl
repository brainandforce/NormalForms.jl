"""
    NormalForms.eye(T::Type, sz)

Generates an identity matrix of size `sz` with elements of type `T`.
"""
eye(T::Type, sz) = Matrix{T}(LinearAlgebra.I(sz))

"""
    NormalForms.eye(M::AbstractMatrix, dim)

Generates an identity matrix of the same dimension as dimension `dim` of `M`.
"""
eye(M::AbstractMatrix, dim) = typeof(M)(collect(UniformScaling(one(eltype(M)))(size(M, dim))))
eye(M::Diagonal, dim = 1) = Diagonal(ones(eltype(M), size(M,dim)))
# Collect is needed for Julia 1.6 because Adjoint(::Diagonal) is undefined

"""
    NormalForms.detb!(M::AbstractMatrix{T<:Number}) -> T

Calculates the determinant of a matrix using the [Bareiss 
algorithm](https://en.wikipedia.org/wiki/Bareiss_algorithm) using in-place operations. This
algorithm returns exact integer determinants for integer matrices, which is useful when checking
whether a matrix is unimodular.

This is a verbatim reimplementation of `LinearAlgebra.det_bareiss()` used to ensure compatibility
with Julia 1.6, which does not have this function.
```
"""
function detb!(M::AbstractMatrix)
    exactdiv(a,b) = a/b
    exactdiv(a::Integer, b::Integer) = div(a,b)
    n = LinearAlgebra.checksquare(M)
    sign, prev = Int8(1), one(eltype(M))
    for i in 1:n-1
        if iszero(M[i,i]) # swap with another col to make nonzero
            swapto = findfirst(!iszero, @view M[i,i+1:end])
            isnothing(swapto) && return zero(prev)
            sign = -sign
            Base.swapcols!(M, i, i + swapto)
        end
        for k in i+1:n, j in i+1:n
            M[j,k] = exactdiv(M[j,k]*M[i,i] - M[j,i]*M[i,k], prev)
        end
        prev = M[i,i]
    end
    return sign * M[end,end]
end

"""
    NormalForms.detb(M::AbstractMatrix) -> eltype(M)

Calculates the determinant of a matrix using the [Bareiss
Algorithm](https://en.wikipedia.org/wiki/Bareiss_algorithm) without modifying `M`.
"""
detb(M::AbstractMatrix) = detb!(copy(M))

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
    NormalForms.gcd_kb(x::Integer, y::Integer) -> NTuple{3,promote_typeof(x,y)}

Calculates the gcd `r` and Bézout coefficients `(p, q)` with the Kannan-Bachem modification. This
ensures that `q ≤ max(x, y)`.

If `x == y`, this function returns `(x, 1, 0)` - this differs from Julia's `gcdx()`, which returns
`(x, 0, 1)`. This modification is needed to prevent `snf_ma!()` from looping infinitely when
zeroing rows and columns.
"""
function gcd_kb(x::Integer, y::Integer)
    x == y && return (x, 1, 0)
    (r,p,q) = gcdx(x,y)
    (b,a) = minmax(x,y)
    if abs(q) > abs(a)
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
@inline function find_pivot(A::AbstractMatrix, k::Integer)
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
    NormalForms.is_row_zero_after(A::AbstractMatrix, k::Integer) -> Bool
    NormalForms.is_column_zero_after(A::AbstractMatrix, k::Integer) -> Bool

Checks if all entries of row or column `k` superseding `A[k,k]` are zero.
"""
is_row_zero_after(A::AbstractMatrix, k::Integer) = all(iszero(A[k,k+1:end]))
is_col_zero_after(A::AbstractMatrix, k::Integer) = all(iszero(A[k+1:end,k]))
@doc (@doc is_row_zero_after) is_col_zero_after

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
    A::AbstractMatrix,
    U::AbstractMatrix,
    k::Integer,
    ki::Integer = k
)
    # This can safely be skipped if k is the last element
    k < last(axes(A,2)) && for j in axes(A,2)[k+1:end]
        # Take the GCD with the Kannen-Bachem modificaion
        (d, p, q) = gcd_kb(A[ki, k], A[ki, j])
        Akkd, Akjd = div(A[ki, k], d), div(A[ki, j], d)
        for i in axes(A,1)
            # Make copies for later reference (otherwise the algorithm breaks)
            Aik, Aij = A[i,k], A[i,j]
            A[i,k] =  Aik * p + Aij * q
            A[i,j] = -Aik * Akjd + Aij * Akkd
        end
        # Modify the unimodular matrix correspondingly
        for i in axes(U,1)
            Uik, Uij = U[i,k], U[i,j]
            U[i,k] =  Uik * p + Uij * q
            U[i,j] = -Uik * Akjd + Uij * Akkd
        end
        #= This is faster than the otherwise equivalent code below:
        Ak, Aj = A[:,k], A[:,j]
        @views A[:,k] =  Ak * p + Aj * q
        @views A[:,j] = -Ak * Akjd + Aj * Akkd
        Uk, Uj = U[:,k], U[:,j]
        @views U[:,k] =  Uk * p + Uj * q
        @views U[:,j] = -Uk * Akjd + Uj * Akkd
        (It was also faster than using a matrix multiplication)
        @assert isunimodular(U) =#
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
    A::AbstractMatrix,
    U::AbstractMatrix,
    k::Integer,
    ki::Integer = k
)
    # This can safely be skipped if k is the last element
    k < last(axes(A,1)) && for j in axes(A,1)[k+1:end]
        (d, p, q) = gcd_kb(A[k, ki], A[j, ki])
        Akkd, Ajkd = div(A[k, ki], d), div(A[j, ki], d)
        for i in axes(A,2)
            Aki, Aji = A[k,i], A[j,i]
            A[k,i] =  Aki * p + Aji * q
            A[j,i] = -Aki * Ajkd + Aji * Akkd
        end
        for i in axes(U,2)
            Uki, Uji = U[k,i], U[j,i]
            U[k,i] =  Uki * p + Uji * q
            U[j,i] = -Uki * Ajkd + Uji * Akkd
        end
        #= This is faster than the otherwise equivalent code below:
        Ak, Aj = A[k,:], A[j,:]
        @views A[k,:] =  Ak * p + Aj * q
        @views A[j,:] = -Ak * Ajkd + Aj * Akkd
        Uk, Uj = U[k,:], U[j,:]
        @views U[k,:] =  Uk * p + Uj * q
        @views U[j,:] = -Uk * Ajkd + Uj * Akkd
        (It was also faster than using a matrix multiplication) =#
    end
end

"""
    NormalForms.reduce_cols_off_diagonal!(
        A::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        R::RoundingMode
    )

Performs an off-diagonal reduction of elements in `A` in-place by columns. This also ensures that
diagonal elements are positive by multiplying the columns by -1 if needed. Changes to `A` are
tracked in the matrix `U`, which will only undergo unimodular transforms.

If the `RoundingMode` is set to `RoundUp` or `NegativeOffDiagonal` (the default), the off-diagonal
elements will be negative. If set to `RoundDown` or `PositiveOffDiagonal`, they will be positive.
"""
@inline function reduce_cols_off_diagonal!(
    A::AbstractMatrix,
    U::AbstractMatrix,
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

"""
    NormalForms.zero_row_and_col!(
        A::AbstractMatrix,
        U::AbstractMatrix,
        V::AbstractMatrix,
        k::Integer,
        ki::Integer = k
    )

Zeros row and column `k` of matrix `A`, tracking the changes in left unimodular matrix `U` and right
unimodular matrix `V`.
"""
function zero_row_and_col!(
    A::AbstractMatrix,
    U::AbstractMatrix,
    V::AbstractMatrix,
    k::Integer,
    ki::Integer = k
)
    # The cycle:
    # Zero the row, reduce the column, zero the column, reduce the row
    zero_row!(A, V, k, ki)
    @assert is_row_zero_after(A, k)
    for n in axes(A, 1)[k+1:end]
        mul = div(A[n,k], A[k,k], RoundToZero)
        A[n,:] .-= mul .* A[k,:]
        U[n,:] .-= mul .* U[k,:]
    end
    zero_col!(A, U, k, ki)
    for n in axes(A, 2)[k+1:end]
        mul = div(A[k,n], A[k,k], RoundToZero)
        A[:,n] .-= mul .* A[:,k]
        V[:,n] .-= mul .* V[:,k]
    end
    @assert is_col_zero_after(A, k)
    @assert is_row_zero_after(A, k)
end

"""
    NormalForms.enforce_divisibility!(
        A::AbstractMatrix,
        U::AbstractMatrix,
        V::AbstractMatrix,
        k::Integer
    )

Ensures that diagonal element `k` is a multiple of diagonal element `k-1` by mutating `A` as well
as the pre-multiplied matrix `U` and post-multiplied matrix `V`.
"""
function enforce_divisibility!(A::AbstractMatrix, U::AbstractMatrix, V::AbstractMatrix, k::Integer)
    # No need to act if there are no previous diagonal entries
    all(isequal(k), first.(axes(A))) && return nothing
    # We want to make A[k-1, k-1] == d
    # ORDER MATTERS HERE: A[k,k] MUST BE SECOND
    (d,p,q) = gcdx(A[k-1,k-1], A[k,k])
    # TODO: does this need to be a while loop?
    while A[k,k] % A[k-1,k-1] != 0
        # Swap the columns to make the diagonal entries zero
        swapcols!(A, k-1, k)
        swapcols!(V, k-1, k)
        # Use the Bezout coefficients to make A[k-1,k-1] == d
        for i in axes(A,2)
            A[k-1,i] += q * A[k,i]
            V[i,k-1] += p * V[i,k]
        end
        for i in axes(A,1)
            A[i,k-1] += p * A[i,k]
            U[k-1,i] += q * U[k,i]
        end
        #= This is faster than the otherwise equivalent code below
        @views A[k-1,:] += q * A[k,:]
        @views U[k-1,:] += q * U[k,:]
        @views A[:,k-1] += p * A[:,k]
        @views V[:,k-1] += p * V[:,k]
        =#
        # Zero the off-diagonal elements
        zero_col!(A, U, k-1)
        zero_row!(A, V, k-1)
    end
    # Make any off-diagonal elmements positive if they went negative
    # The branchless approach here has the same efficiency as a branched one
    for n in (-1:0) .+ k
        sg = sign(A[n,n])
        for i in axes(A,1)
            U[n,i] *= sg
        end
        A[n,n] *= sg
    end
end

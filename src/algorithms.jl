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
    if q > x && q > y
        (a,b) = ((x,y),) |> (minimum, maximum)
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
function find_pivot(A::AbstractMatrix{<:Integer}, k::Integer)
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

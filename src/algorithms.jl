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

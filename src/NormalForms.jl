module NormalForms

using LinearAlgebra
using Requires
import LinearAlgebra: det_bareiss as detb
import Base: promote_typeof

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

export isunimodular

include("hermite.jl")
export AbstractHermite, RowHermite, ColumnHermite
export NegativeOffDiagonal, PositiveOffDiagonal
export hnfc!, hnfc, hnfr!, hnfr
include("smith.jl")
export Smith
export snf!, snf

function __init__()
    # Contains methods specific to SMatrix, which can't be mutated in-place
    @require StaticArrays="90137ffa-7385-5640-81b9-e52037218182" begin
        include("staticarrays.jl")
    end
end

end

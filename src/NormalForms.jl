module NormalForms

using LinearAlgebra
using Requires
import Base.setindex

# Efficient algorithm for calculating determinants of integer matrices
const detb = LinearAlgebra.det_bareiss

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

export isunimodular

include("hermite.jl")
export AbstractHermite, RowHermite, ColumnHermite
export NegativeOffDiagonal, PositiveOffDiagonal
export hnfc!, hnfc, hnfr!, hnfr
include("smith.jl")
export Smith

function __init__()
    # Contains methods specific to SMatrix, which can't be mutated in-place
    @require StaticArrays="90137ffa-7385-5640-81b9-e52037218182" begin
        include("staticarrays.jl")
    end
end

end

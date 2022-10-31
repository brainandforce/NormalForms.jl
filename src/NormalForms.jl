module NormalForms

using LinearAlgebra

# Efficient algorithm for calculating determinants of integer matrices
const detb = LinearAlgebra.det_bareiss

"""
    isunimodular(M::AbstractMatrix)

Returns `true` if a matrix is unimodular. A unimodular matrix is a square integer matrix with A
determinant equal to 1 or -1.
"""
isunimodular(M::AbstractMatrix{<:Integer}) = issquare(M) && isone(abs(detb(M)))

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
include("smith.jl")
export Smith

end

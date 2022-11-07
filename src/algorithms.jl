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

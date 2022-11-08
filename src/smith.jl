"""
    Smith{T,M<:AbstractMatrix{T}} <: Factorization{T}

Describes the result of the decomposition of an integer matrix `A` into its Smith normal form `S`,
which is diagonal, and unimodular matrices `U` and `V` such that `S == U*A*V`, or equivalently,
`A == U\\S/V`.
"""
struct Smith{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}
    S::M
    U::M
    V::M
    info::LinearAlgebra.BlasInt
    function Smith(
        S::AbstractMatrix{<:Integer},
        U::AbstractMatrix{<:Integer},
        V::AbstractMatrix{<:Integer},
        info::Integer = 0
    )
        @assert isunimodular(U) "U is not unimodular."
        @assert isunimodular(V) "V is not unimodular."
        @assert size(S,1) == size(U,2) string(
            "S and U have incompatible dimensions: ",
            join(size(S), '×'), " for S and ", join(size(U), '×'), " for U"
        )
        @assert size(S,2) == size(V,1) string(
            "S and V have incompatible dimensions: ",
            join(size(S), '×'), " for S and ", join(size(V), '×'), " for V"
        )
        M = promote_typeof(S, U, V)
        return new{eltype(M),M}(S, U, V, info)
    end
end

# Component destructuring
Base.iterate(F::Smith) = (F.S, Val{:U}())
Base.iterate(F::Smith, ::Val{:S}) = iterate(F)
Base.iterate(F::Smith, ::Val{:U}) = (F.U, Val{:V}())
Base.iterate(F::Smith, ::Val{:V}) = (F.V, nothing)
Base.iterate(F::Smith, ::Nothing) = nothing

issuccess(F::Smith) = iszero(F.info)

function Base.summary(io::IO, F::Smith)
    print(io, join(size(F.S), '×'), ' ', typeof(F), ":")
end

function Base.show(io::IO, mime::MIME"text/plain", F::Smith)
    if issuccess(F)
        summary(io, F)
        println(io, "\nSmith normal form:")
        Base.show(io, mime, F.S)
        println(io, "\nLeft unimodular factor:")
        Base.show(io, mime, F.U)
        println(io, "\nRight unimodular factor:")
        Base.show(io, mime, F.V)
    else
        print(io, "Failed factorization of type $(typeof(F))")
    end
end

function snf_ma!(A::AbstractMatrix{<:Integer})
    # Create the left and right unimodular matrices
    U = eye(A,1)
    V = eye(A,2)
    # Convert to a transpose or adjoint if needed
    if A isa Adjoint
        U = U'
        V = V'
    elseif A isa Transpose
        U = transpose(U)
        V = transpose(V)
    end
    # Loop through the smallest dimension of A
    for k in minimum(axes(A))
        pivot = find_pivot(A,k)
        # It's not entirely clear to me what ki does.
        if iszero(pivot)
            # Set ki to the first nonzero column member
            ki = findfirst(!iszero, @view A[k+1:end, k])
            # If all are zero, break from the row loop (matrix is truly rank deficient)
            # Otherwise, add the found value to k to compensate for the view's indices
            isnothing(ki) ? continue : (ki += k) # what does it mean when ki > k...
        else
            ki = k
            swapcols!(A, pivot, k)
            swapcols!(V, pivot, k)
        end
        # Loop until all off-diagonal elements are zero
        counter = 0
        while !all(iszero, (A[k,k+1:end], A[k+1:end,k]))
            # Zero off-diagonal elements across the row, then the column
            zero_row!(A, V, k, ki)
            zero_col!(A, U, k, ki)
            @info "A = " * repr("text/plain", A) * 
                "\nA[$k,$(k+1):end] = " * repr(A[k,k+1:end]) *
                "\nA[$(k+1):end,$k] = " * repr(A[k+1:end,k])
            counter += 1
            @info "On iteration $counter (k = $k)"
            if counter == 32
                error("We got stuck.")
            end
        end
        # Perform some new operations to make A[k,k] divisible by A[k-1,k-1]
        enforce_divisibility!(A, U, V, k)
    end
    return (A, U, V, 0)
end

"""
    snf(M::AbstractMatrix{<:Integer}) -> Smith{eltype(M),M}

Calculates the Smith normal form of an integer matrix in-place.
"""
snf!(M::AbstractMatrix{<:Integer}) = Smith(snf_ma!(M)...)

"""
    snf(M::AbstractMatrix{<:Integer}) -> Smith{eltype(M),M}

Calculates the Smith normal form of an integer matrix, returning a copy of the original matrix.
"""
snf(M::AbstractMatrix{<:Integer}) = snf!(deepcopy(M))

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
        println(io, "\nHermite normal form:")
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
        # Set k as a default pivot
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
        # It's not entirely clear to me what ki does.
        if iszero(pivot)
            # Set ki to the first nonzero column member
            ki = findfirst(!iszero, @view A[k+1:end, k])
            # If all are zero, break from the row loop (matrix is truly rank deficient)
            # Otherwise, add the found value to k to compensate for the view's indices
            isnothing(ki) ? continue : (ki += k) # what does it mean when ki > k...
        else
            ki = k
            # Permute rows according to the pivot if needed
            if pivot != k
                # TODO: Is this the right set of operations?
                A[:,[pivot,k]] = A[:,[k,pivot]]
                # We're only pivoting by column: U[[pivot,k],:] = U[[k,pivot],:]
                V[:,[pivot,k]] = V[:,[k,pivot]]
            end
        end
        # Loop until all off-diagonal elements are zero
        counter = 0
        while !all(iszero, (A[k,k+1:end], A[k+1:end,k]))
            # Zero the off-diagonal elements across the row: A[k,k+1:end]
            # @info "row = $row, col = $col"
            k < size(A,2) && for j in axes(A,2)[k+1:end]
                # Generate the matrix elements for this transform
                Akk, Akj = A[ki, k], A[ki, j]
                (d, p, q) = gcdx(Akk, Akj)
                # Rounding is irrelevant since d is the gcd
                Akkd, Akjd = div(Akk, d), div(Akj, d)
                # Mutating A to zero the upper off-diagonal elements
                Ak, Aj = A[:,k], A[:,j]
                A[:,k] =  Ak * p + Aj * q
                A[:,j] = -Ak * Akjd + Aj * Akkd
                # Mutating V as the matching unimodular matrix
                Vk, Vj = V[:,k], V[:,j]
                V[:,k] =  Vk * p + Vj * q
                V[:,j] = -Vk * Akjd + Vj * Akkd
                # @assert isunimodular(V)
            end
            #=
            @info "A = " * repr("text/plain", A) * 
                "\nA[$k,$(k+1):end] = " * repr(A[k,k+1:end]) *
                "\nA[$(k+1):end,$k] = " * repr(A[k+1:end,k])
            =#
            # Now zero them across the column: A[k+1:end,k]
            k < size(A,1) && for j in axes(A,1)[k+1:end]
                # Generate the matrix elements for this transform
                Akk, Ajk = A[k, ki], A[j, ki]
                (d, p, q) = gcdx(Akk, Ajk)
                # Rounding is irrelevant since d is the gcd
                Akkd, Ajkd = div(Akk, d), div(Ajk, d)
                # @info string("D_row = ", repr("text/plain", [p q; Ajkd Akkd]))
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
            @info "A = " * repr("text/plain", A) * 
                "\nA[$k,$(k+1):end] = " * repr(A[k,k+1:end]) *
                "\nA[$(k+1):end,$k] = " * repr(A[k+1:end,k])
            counter += 1
            @info "On iteration $counter (k = $k)"
            if counter == 32
                error("We got stuck.")
            end
        end
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

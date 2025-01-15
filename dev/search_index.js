var documenterSearchIndex = {"docs":
[{"location":"notes/#Package-notes-and-FAQ","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"","category":"section"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"Some of the content here may seem obvious, but these are pitfalls that were discovered during package development that led to some design choices and other unusual findings that are rationalized or explained here.","category":"page"},{"location":"notes/#Return-type-of-Transpose-and-Adjoint-matrices","page":"Package notes and FAQ","title":"Return type of Transpose and Adjoint matrices","text":"","category":"section"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"When a Smith or Hermite normal form calculation is performed on a Transpose or Adjoint, the unimodular factors are of the same type. The transpose() and adjoint() functions can be applied to a factorization to obtain a transposed representation.","category":"page"},{"location":"notes/#Inequivalence-of-transposition-order-with-Smith-normal-form-calculations","page":"Package notes and FAQ","title":"Inequivalence of transposition order with Smith normal form calculations","text":"","category":"section"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"In general, hnfc(M') == hnfr(M)', and hnfr(M') == hnfc(M)'. It is also true that  snf(M).S == snf(M')'.S == snf(M)'.S. However, snf(M').U is not guaranteed to be equal to snf(M)'.U, and the same goes for snf(M').V and snf(M)'.V! This is because the choice of unimodular factors is not unique with respect to the rounding mode of the Euclidean division steps.","category":"page"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"What is guaranteed is that snf(M).S == snf(M').S == snf(M)'.S, snf(M)'.S == snf(M').S == snf(M').V * M * snf(M').U == snf(M)'.V * M' * snf(M)'.U. In the future, we may modify the algorithm to guarantee snf(M') == snf(M)'.","category":"page"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"note: Note\nWe've generally used the adjoint (prime) operator, but everything stated for adjoints here also holds for transposes.","category":"page"},{"location":"notes/#Return-type-of-Diagonal-matrices","page":"Package notes and FAQ","title":"Return type of Diagonal matrices","text":"","category":"section"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"You might notice that if you perform a Smith normal form calculation on a Diagonal matrix, the result is a Smith{T, Matrix{T}}, instead of Smith{T, Diagonal{T, Vector{T}}}. However, this is not the case for Hermite normal form calculations, which will return an AbstractHermite{T, Diagonal{T, Vector{T}}}:","category":"page"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"julia> M = Diagonal([7,12,6])\n3×3 Diagonal{Int64, Vector{Int64}}:\n 7   ⋅  ⋅\n ⋅  12  ⋅\n ⋅   ⋅  6\n\njulia> snf(M)\n3×3 Smith{Int64, Matrix{Int64}}:\nSmith normal form:\n3×3 Matrix{Int64}:\n 1  0   0\n 0  6   0\n 0  0  84\nLeft unimodular factor:\n3×3 Matrix{Int64}:\n  1   3  0\n 12  35  1\n 12  35  0\nRight unimodular factor:\n3×3 Matrix{Int64}:\n -5  0   36\n  1  0   -7\n  0  1  -14\n\njulia> hnfc(M)\n3×3 ColumnHermite{Int64, Diagonal{Int64, Vector{Int64}}}:\nHermite normal form:\n3×3 Diagonal{Int64, Vector{Int64}}:\n 7   ⋅  ⋅\n ⋅  12  ⋅\n ⋅   ⋅  6\nUnimodular factor:\n3×3 Diagonal{Int64, Vector{Int64}}:\n 1  ⋅  ⋅\n ⋅  1  ⋅\n ⋅  ⋅  1\n\njulia> hnfr(M)\n3×3 RowHermite{Int64, Diagonal{Int64, Vector{Int64}}}:\nHermite normal form:\n3×3 Diagonal{Int64, Vector{Int64}}:\n 7   ⋅  ⋅\n ⋅  12  ⋅\n ⋅   ⋅  6\nUnimodular factor:\n3×3 Diagonal{Int64, Vector{Int64}}:\n 1  ⋅  ⋅\n ⋅  1  ⋅\n ⋅  ⋅  1\n","category":"page"},{"location":"notes/","page":"Package notes and FAQ","title":"Package notes and FAQ","text":"A diagonal matrix automatically satisfies all of the requirements for being in Hermite normal form, so the unimodular factor is just an identity matrix, which can be stored in a Diagonal representation. However, it does not automatically satisfy the divisibility requirements for a  matrix in Smith normal form. Therefore, the matrix is converted to a dense matrix so that the unimodular factors, which themselves must be stored as dense matrices, can be permuted in the usual way.","category":"page"},{"location":"api/internal/#Internal-APIs","page":"Internal","title":"Internal APIs","text":"","category":"section"},{"location":"api/internal/","page":"Internal","title":"Internal","text":"CurrentModule = NormalForms","category":"page"},{"location":"api/internal/","page":"Internal","title":"Internal","text":"Most users should not need to read this section: it's primarily intended for those who would like to contribute to the package.","category":"page"},{"location":"api/internal/","page":"Internal","title":"Internal","text":"Docstrings for internal functions and other objects begin with NormalForms, in contrast with the public API.","category":"page"},{"location":"api/internal/#Factorizations","page":"Internal","title":"Factorizations","text":"","category":"section"},{"location":"api/internal/","page":"Internal","title":"Internal","text":"NormalForms.hnf_ma!\nNormalForms.snf_ma!","category":"page"},{"location":"api/internal/#NormalForms.hnf_ma!","page":"Internal","title":"NormalForms.hnf_ma!","text":"NormalForms.hnf_ma!(\n    A::AbstractMatrix{<:Integer},\n    R::RoundingMode = NegativeOffDiagonal\n) -> Tuple{typeof(A), typeof(A), LinearAlgebra.BlasInt}\n\nA modified version of HermiteNormalForm._hnf_like!() that calculates the column-style Hermite normal form of a matrix in-place, returning A, the unimodular matrix U such that H == A*U, and the factorization info (zero if success - matching LinearAlgebra factorizations).\n\nThe original implementation of HermiteNormalForm._hnf_like!() was written by Yingbo Ma, and may be found here: https://github.com/YingboMa/HermiteNormalForm.jl/\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.snf_ma!","page":"Internal","title":"NormalForms.snf_ma!","text":"NormalForms.snf_ma!(A::AbstractMatrix{<:Integer})\n    -> Tuple{typeof(A), typeof(A), typeof(A), LinearAlgebra.BlasInt}\n\nA modified version of HermiteNormalForm._hnf_like!() that calculates the Smith normal form of a matrix in-place, returning A, the unimodular matrices U and V such that S == U*A*V, and the factorization info (zero if success - matching LinearAlgebra factorizations).\n\nThe original implementation of HermiteNormalForm._hnf_like!() was written by Yingbo Ma, and may be found here: https://github.com/YingboMa/HermiteNormalForm.jl/\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#Algorithms","page":"Internal","title":"Algorithms","text":"","category":"section"},{"location":"api/internal/","page":"Internal","title":"Internal","text":"NormalForms.eye\nNormalForms.gcd_kb\nNormalForms.detb!\nNormalForms.detb\nNormalForms.find_pivot\nNormalForms.is_row_zero_after\nNormalForms.is_col_zero_after\nNormalForms.zero_row!\nNormalForms.zero_col!\nNormalForms.reduce_cols_off_diagonal!\nNormalForms.zero_row_and_col!\nNormalForms.enforce_divisibility!","category":"page"},{"location":"api/internal/#NormalForms.eye","page":"Internal","title":"NormalForms.eye","text":"NormalForms.eye(T::Type, sz)\n\nGenerates an identity matrix of size sz with elements of type T.\n\n\n\n\n\nNormalForms.eye(M::AbstractMatrix, dim)\n\nGenerates an identity matrix of the same dimension as dimension dim of M.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.gcd_kb","page":"Internal","title":"NormalForms.gcd_kb","text":"NormalForms.gcd_kb(x::Integer, y::Integer) -> NTuple{3,promote_typeof(x,y)}\n\nCalculates the gcd r and Bézout coefficients (p, q) with the Kannan-Bachem modification. This ensures that abs(q) ≤ max(abs(x), abs(y)).\n\nIf x == y, this function returns (x, 1, 0) - this differs from Julia's gcdx(), which returns (x, 0, 1). This modification is needed to prevent snf_ma!() from looping infinitely when zeroing rows and columns.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.detb!","page":"Internal","title":"NormalForms.detb!","text":"NormalForms.detb!(M::AbstractMatrix{T<:Number}) -> T\n\nCalculates the determinant of a matrix using the Bareiss  algorithm using in-place operations. This algorithm returns exact integer determinants for integer matrices, which is useful when checking whether a matrix is unimodular.\n\nThis is a verbatim reimplementation of LinearAlgebra.det_bareiss() used to ensure compatibility with Julia 1.6, which does not have this function. ```\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.detb","page":"Internal","title":"NormalForms.detb","text":"NormalForms.detb(M::AbstractMatrix) -> eltype(M)\n\nCalculates the determinant of a matrix using the Bareiss Algorithm without modifying M.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.find_pivot","page":"Internal","title":"NormalForms.find_pivot","text":"NormalForms.find_pivot(A::AbstractMatrix{<:Integer}, k::Integer)\n\nFinds a suitable pivot element for the Smith or Hermite normal form algorithms. If none can be found, it returns 0.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.is_row_zero_after","page":"Internal","title":"NormalForms.is_row_zero_after","text":"NormalForms.is_row_zero_after(A::AbstractMatrix, k::Integer) -> Bool\nNormalForms.is_column_zero_after(A::AbstractMatrix, k::Integer) -> Bool\n\nChecks if all entries of row or column k superseding A[k,k] are zero.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.is_col_zero_after","page":"Internal","title":"NormalForms.is_col_zero_after","text":"NormalForms.is_row_zero_after(A::AbstractMatrix, k::Integer) -> Bool\nNormalForms.is_column_zero_after(A::AbstractMatrix, k::Integer) -> Bool\n\nChecks if all entries of row or column k superseding A[k,k] are zero.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.zero_row!","page":"Internal","title":"NormalForms.zero_row!","text":"NormalForms.zero_row!(\n    A::AbstractMatrix{<:Integer},\n    U::AbstractMatrix{<:Integer},\n    k::Integer\n    ki::Integer = k\n)\n\nZeroes the elements along the row A[k,k+1:end] (all elements right of A[k,k]). Changes to A are tracked in the matrix U, which will only undergo unimodular transforms.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.zero_col!","page":"Internal","title":"NormalForms.zero_col!","text":"NormalForms.zero_col!(\n    A::AbstractMatrix{<:Integer},\n    U::AbstractMatrix{<:Integer},\n    k::Integer\n    ki::Integer = k\n)\n\nZeroes the elements along the column A[k+1:end,k] (all elements below A[k,k]). Changes to A are tracked in the matrix U, which will only undergo unimodular transforms.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.reduce_cols_off_diagonal!","page":"Internal","title":"NormalForms.reduce_cols_off_diagonal!","text":"NormalForms.reduce_cols_off_diagonal!(\n    A::AbstractMatrix{<:Integer},\n    U::AbstractMatrix{<:Integer},\n    R::RoundingMode\n)\n\nPerforms an off-diagonal reduction of elements in A in-place by columns. This also ensures that diagonal elements are positive by multiplying the columns by -1 if needed. Changes to A are tracked in the matrix U, which will only undergo unimodular transforms.\n\nIf the RoundingMode is set to RoundUp or NegativeOffDiagonal (the default), the off-diagonal elements will be negative. If set to RoundDown or PositiveOffDiagonal, they will be positive.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.zero_row_and_col!","page":"Internal","title":"NormalForms.zero_row_and_col!","text":"NormalForms.zero_row_and_col!(\n    A::AbstractMatrix,\n    U::AbstractMatrix,\n    V::AbstractMatrix,\n    k::Integer,\n    ki::Integer = k\n)\n\nZeros row and column k of matrix A, tracking the changes in left unimodular matrix U and right unimodular matrix V.\n\n\n\n\n\n","category":"function"},{"location":"api/internal/#NormalForms.enforce_divisibility!","page":"Internal","title":"NormalForms.enforce_divisibility!","text":"NormalForms.enforce_divisibility!(\n    A::AbstractMatrix,\n    U::AbstractMatrix,\n    V::AbstractMatrix,\n    k::Integer\n)\n\nEnsures that diagonal element k is a multiple of diagonal element k-1 by mutating A as well as the pre-multiplied matrix U and post-multiplied matrix V.\n\n\n\n\n\n","category":"function"},{"location":"#NormalForms","page":"Home","title":"NormalForms","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides functions for calculating the Hermite normal form (hnfr!(), hnfr(),  hnfc!(), hnfc()), and the Smith normal form (snf!(), snf!()), as well as factorization types RowHermite, ColumnHermite, and Smith, which are subtypes of LinearAlgebra.Factorization.","category":"page"},{"location":"api/public/#Public-API","page":"Public","title":"Public API","text":"","category":"section"},{"location":"api/public/","page":"Public","title":"Public","text":"CurrentModule = NormalForms","category":"page"},{"location":"api/public/","page":"Public","title":"Public","text":"This consists of all of the exported functions and objects available to the end user.","category":"page"},{"location":"api/public/#Factorizations","page":"Public","title":"Factorizations","text":"","category":"section"},{"location":"api/public/","page":"Public","title":"Public","text":"hnfc!\nhnfc\nhnfr!\nhnfr\nsnf!\nsnf","category":"page"},{"location":"api/public/#NormalForms.hnfc!","page":"Public","title":"NormalForms.hnfc!","text":"hnfc!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)\n    -> ColumnHermite{eltype(M),typeof(M)}\n\nCalculates the column Hermite normal form of the integer matrix M in-place, returning a ColumnHermite describing the elements of the factorization.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#NormalForms.hnfc","page":"Public","title":"NormalForms.hnfc","text":"hnfc(M::AbstractMatrix{T<:Integer}) -> ColumnHermite{T,typeof(M)}\n    -> ColumnHermite{eltype(M),typeof(M)}\n\nCalculates the column Hermite normal form of the integer matrix M in-place, returning a ColumnHermite describing the elements of the factorization. Unlike hnfc!(), this function creates a copy of the input matrix rather than modifying it in-place.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#NormalForms.hnfr!","page":"Public","title":"NormalForms.hnfr!","text":"hnfr!(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)\n    -> RowHermite{eltype(M),typeof(M)}\n\nCalculates the row Hermite normal form of the integer matrix M in-place, returning a RowHermite describing the elements of the factorization.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#NormalForms.hnfr","page":"Public","title":"NormalForms.hnfr","text":"hnfr(M::AbstractMatrix{<:Integer}, R::RoundingMode = NegativeOffDiagonal)\n    -> RowHermite{eltype(M),typeof(M)}\n\nCalculates the row Hermite normal form of the integer matrix M in-place, returning a RowHermite describing the elements of the factorization. Unlike hnfr!(), this function creates a copy of the input matrix rather than modifying it in-place.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#NormalForms.snf!","page":"Public","title":"NormalForms.snf!","text":"snf(M::AbstractMatrix{<:Integer}) -> Smith{eltype(M),M}\n\nCalculates the Smith normal form of an integer matrix in-place.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#NormalForms.snf","page":"Public","title":"NormalForms.snf","text":"snf(M::AbstractMatrix{<:Integer}) -> Smith{eltype(M),M}\n\nCalculates the Smith normal form of an integer matrix, returning a copy of the original matrix.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#Helper-functions","page":"Public","title":"Helper functions","text":"","category":"section"},{"location":"api/public/","page":"Public","title":"Public","text":"isunimodular","category":"page"},{"location":"api/public/#NormalForms.isunimodular","page":"Public","title":"NormalForms.isunimodular","text":"isunimodular(M::AbstractMatrix)\n\nReturns true if a matrix is unimodular. A unimodular matrix is a square integer matrix with A determinant equal to 1 or -1.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#Factorization-objects","page":"Public","title":"Factorization objects","text":"","category":"section"},{"location":"api/public/","page":"Public","title":"Public","text":"AbstractHermite\nRowHermite\nColumnHermite\nSmith","category":"page"},{"location":"api/public/#NormalForms.AbstractHermite","page":"Public","title":"NormalForms.AbstractHermite","text":"AbstractHermite{T<:Integer,M<:AbstractMatrix{T}} <: Factorization{T}\n\nSupertype for the result of Hermite normal form calculations.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#NormalForms.RowHermite","page":"Public","title":"NormalForms.RowHermite","text":"RowHermite{T<:Integer,M<:AbstractMatrix{T}}\n\nDescribes the result of a decomposition of an integer matrix A into its row-style Hermite normal form H, which is upper triangular, and a unimodular matrix U, such that H == U*A, or equivalently, A = U\\H.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#NormalForms.ColumnHermite","page":"Public","title":"NormalForms.ColumnHermite","text":"ColumnHermite{T<:Integer,M<:AbstractMatrix{T}}\n\nDescribes the result of a decomposition of an integer matrix A into its row-style Hermite normal form H, which is lower triangular, and a unimodular matrix U, such that H == A*U, or equivalently, A == H/U.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#NormalForms.Smith","page":"Public","title":"NormalForms.Smith","text":"Smith{T,M<:AbstractMatrix{T}} <: Factorization{T}\n\nDescribes the result of the decomposition of an integer matrix A into its Smith normal form S, which is diagonal, and unimodular matrices U and V such that S == U*A*V, or equivalently, A == U\\S/V.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#Reduction-modes","page":"Public","title":"Reduction modes","text":"","category":"section"},{"location":"api/public/","page":"Public","title":"Public","text":"NegativeOffDiagonal\nPositiveOffDiagonal","category":"page"},{"location":"api/public/#NormalForms.NegativeOffDiagonal","page":"Public","title":"NormalForms.NegativeOffDiagonal","text":"NegativeOffDiagonal\n\nAlias for RoundUp. In NormalForms.reduce_cols_off_diagonal!(), this makes all off-diagonal elements the smallest negative value possible.\n\n\n\n\n\n","category":"constant"},{"location":"api/public/#NormalForms.PositiveOffDiagonal","page":"Public","title":"NormalForms.PositiveOffDiagonal","text":"PositiveOffDiagonal\n\nAlias for RoundDown. In NormalForms.reduce_cols_off_diagonal!(), this makes all off-diagonal elements the smallest positive value possible.\n\n\n\n\n\n","category":"constant"}]
}

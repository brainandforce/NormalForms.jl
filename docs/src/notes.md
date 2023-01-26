# Package notes and FAQ

Some of the content here may seem obvious, but these are pitfalls that were discovered during
package development that led to some design choices and other unusual findings that are
rationalized or explained here.

## Return type of `Transpose` and `Adjoint` matrices

When a Smith or Hermite normal form calculation is performed on a `Transpose` or `Adjoint`, the
unimodular factors are of the same type. The `transpose()` and `adjoint()` functions can be applied
to a factorization to obtain a transposed representation.

## Inequivalence of transposition order with Smith normal form calculations

In general, `hnfc(M') == hnfr(M)'`, and `hnfr(M') == hnfc(M)'`. It is also true that 
`snf(M).S == snf(M')'.S == snf(M)'.S`. However, `snf(M').U` is not guaranteed to be equal to
`snf(M)'.U`, and the same goes for `snf(M').V` and `snf(M)'.V`! This is because the choice of
unimodular factors is not unique with respect to the rounding mode of the Euclidean division steps.

What is guaranteed is that `snf(M).S == snf(M').S == snf(M)'.S`,
`snf(M)'.S == snf(M').S == snf(M').V * M * snf(M').U == snf(M)'.V * M' * snf(M)'.U`. In the future,
we may modify the algorithm to guarantee `snf(M') == snf(M)'`.

!!! note 
    We've generally used the adjoint (prime) operator, but everything stated for adjoints here
    also holds for transposes.

## Return type of `Diagonal` matrices

You might notice that if you perform a Smith normal form calculation on a `Diagonal` matrix, the
result is a `Smith{T, Matrix{T}}`, instead of `Smith{T, Diagonal{T, Vector{T}}}`. However, this is
not the case for Hermite normal form calculations, which will return an
`AbstractHermite{T, Diagonal{T, Vector{T}}}`:

```julia-repl
julia> M = Diagonal([7,12,6])
3×3 Diagonal{Int64, Vector{Int64}}:
 7   ⋅  ⋅
 ⋅  12  ⋅
 ⋅   ⋅  6

julia> snf(M)
3×3 Smith{Int64, Matrix{Int64}}:
Smith normal form:
3×3 Matrix{Int64}:
 1  0   0
 0  6   0
 0  0  84
Left unimodular factor:
3×3 Matrix{Int64}:
  1   3  0
 12  35  1
 12  35  0
Right unimodular factor:
3×3 Matrix{Int64}:
 -5  0   36
  1  0   -7
  0  1  -14

julia> hnfc(M)
3×3 ColumnHermite{Int64, Diagonal{Int64, Vector{Int64}}}:
Hermite normal form:
3×3 Diagonal{Int64, Vector{Int64}}:
 7   ⋅  ⋅
 ⋅  12  ⋅
 ⋅   ⋅  6
Unimodular factor:
3×3 Diagonal{Int64, Vector{Int64}}:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1

julia> hnfr(M)
3×3 RowHermite{Int64, Diagonal{Int64, Vector{Int64}}}:
Hermite normal form:
3×3 Diagonal{Int64, Vector{Int64}}:
 7   ⋅  ⋅
 ⋅  12  ⋅
 ⋅   ⋅  6
Unimodular factor:
3×3 Diagonal{Int64, Vector{Int64}}:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1

```

A diagonal matrix automatically satisfies all of the requirements for being in Hermite normal form,
so the unimodular factor is just an identity matrix, which can be stored in a `Diagonal`
representation. However, *it does not automatically satisfy the divisibility requirements for a 
matrix in Smith normal form.* Therefore, the matrix is converted to a dense matrix so that the
unimodular factors, which themselves must be stored as dense matrices, can be permuted in the usual
way.

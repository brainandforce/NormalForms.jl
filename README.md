# NormalForms.jl

[![Documentation (stable)][docs-stable-img]][docs-stable-url]
[![Documentation (dev)][docs-dev-img]][docs-dev-url]
[![CI status][ci-status-img]][ci-status-url]
[![Codecov][codecov-img]][codecov-url]
[![Aqua.jl][aqua-img]][aqua-url]

This package allows for the calculation of both the Hermite and Smith normal forms, which are
commonly used throughout crystallography.

Packages to calculate the Smith and Hermite normal forms already exist, but this package provides
several advantages. First, it integrates into pre-existing infrastructure provided by the
`LinearAlgebra` standard library by exporting the `RowHermite`, `ColumnHermite`, and `Smith` types,
which are subtypes of `LinearAlgebra.Factorization`. Second, effort has been made to thoroughly
comment the code so the algorithms can be easily understood. Third, the code is integrated with
StaticArrays.jl, and methods are available for `SMatrix`, which cannot be mutated in-place.

This package only works with Julia 1.7 and above due to the use of `LinearAlgebra.det_bareiss` to
calculate exact integer determinants.

# To-do list

This package is not complete by any means, and the following outstanding issues exist:
  * Smith normal forms are not yet available.
  * Documentation and testing are far from complete.

# See also

[HermiteNormalForm.jl](https://github.com/YingboMa/HermiteNormalForm.jl)

[SmithNormalForm.jl](https://github.com/wildart/SmithNormalForm.jl) (Not available in the General
packge repository)

For those interested in the Hermite and Smith normal forms in the context of abstract algebra:

[Nemo.jl](https://github.com/Nemocas/Nemo.jl)

[Hecke.jl](https://github.com/thofma/Hecke.jl)

[docs-stable-img]:  https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:  https://brainandforce.github.io/NormalForms.jl/stable
[docs-dev-img]:     https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:     https://brainandforce.github.io/NormalForms.jl/dev
[ci-status-img]:    https://github.com/brainandforce/NormalForms.jl/workflows/CI/badge.svg
[ci-status-url]:    https://github.com/brainandforce/NormalForms.jl/actions
[aqua-img]:         https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]:         https://github.com/JuliaTesting/Aqua.jl
[codecov-img]:      https://codecov.io/gh/brainandforce/NormalForms.jl/branch/main/graph/badge.svg
[codecov-url]:      https://codecov.io/gh/brainandforce/NormalForms.jl/

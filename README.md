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

The 0.1 series is compatible with Julia 1.6 and later. All future breaking releases will support the
latest LTS at the time of their release. Therefore, you can expect version 0.2.0 to support Julia
1.10 and later.

# Installation

NormalForms.jl is now in the Julia package registry! To install, just use the package manager:
```
(@v1.6+) pkg> add NormalForms
```
Alternatively: 
```
julia> import Pkg

Pkg.add("NormalForms")
```
You can also add the repo URL, with an optional specifier for the branch you want to track. By
default, it tracks `main`. `release` matches the package state in the General repo.
```
(@v1.6+) pkg> add https://github.com/brainandforce/NormalForms.jl#branchname
```
Alternatively: 
```
julia> import Pkg

Pkg.add("https://github.com/brainandforce/NormalForms.jl", rev="branchname")
```

# To-do list

This package is not in a finished state. The primary block to this is that some matrices cannot be
placed in Hermite or Smith normal form because the unimodular factors fail to be unimodular. As of
now, it's unclear what causes this to happen, but the solution may come in more extensive testing.

# See also

* [HermiteNormalForm.jl](https://github.com/YingboMa/HermiteNormalForm.jl) (This package adapts some
code from here)
* [SmithNormalForm.jl](https://github.com/wildart/SmithNormalForm.jl) (Not available in the General
package repository)

For those interested in the Hermite and Smith normal forms in the context of abstract algebra:

* [Nemo.jl](https://github.com/Nemocas/Nemo.jl)
* [Hecke.jl](https://github.com/thofma/Hecke.jl)

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

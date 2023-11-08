# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## Fixed
  - Issue #13 was not fixed by v0.1.6; it should be fixed now.

## [0.1.6] - 2023-11-05

## Added
  - New methods: `NormalForms.is_row_zero_after`, `NormalForms.is_col_zero_after`, 
`NormalForms.zero_row_and_col!`

### Fixed
  - Row/column zeroing should not loop forever (#13)

## [0.1.5] - 2023-01-25

### Added
  - Factorizations can now be converted to a `Diagonal` representation by calling `Diagonal()`.
  - Factorizations may now be transposed, or converted to an adjoint representation.

### Fixed
  - The method `snf_ma!(::Diagonal)` did not perform any operations, which of course skipped over
the steps that ensure that the diagonal entries meet divisibility requirements. This has been fixed
by converting `Diagonal` matrices to the `Matrix` type before proceeding.

## [0.1.4] - 2023-01-24

### Added
  - All Hermite and Smith normal form calculating functions now accept non-integer inputs, but
convert the input matrix to an integer matrix before proceeding.

### Fixed
  - Diagonal matrices now have their own methods for `hnf_ma!()` and `snf_ma!()`, preventing type
errors and improving efficiency.

### Other
  - Broadened method signature types for many internal algorithms.

## [0.1.3] - 2023-01-20

### Fixed
  - Fixed a bug where matrices with off-diagonal entries equal to diagonal entries may never be
able to become zero off-diagonal.

### Other
  - Updated `NormalForms.eye()` to handle transpose and adjoint matrices more transparently.

## [0.1.2] - 2022-12-12

### Fixed
  - `snf(::SMatrix{D1,D2,T})` now returns its factors as `SMatrix` types.
  - The return type of `snf(::SMatrix{D1,D2,T})` is now `Smith{T,SArray{S,T,2}} where S<:Tuple`.
  - Removed an extraneous `@info` statement from `NormalForms.snf_ma!()`

## [0.1.1] - 2022-11-08

### Added

  - This changelog!
  - Complete documentation of all functions, public and internal
  - Detailed installation instructions in [README.md]
  - Tests for factorization constructors, iteration specification, and some algorithms
  - Bareiss algorithm is now included in package to support Julia releases before 1.7
  - `diag()` definition for all factorizations

### Fixed

  - Broken iterator for `AbstractHermite` types (`RowHermite` and `ColumnHermite`)

### Removed
  
  - Unused `HermiteStyle` dispatch type

## [0.1.0] - 2022-11-07

Initial release - added to Julia package registry

[Unreleased]: https://github.com/brainandforce/NormalForms.jl
[0.1.6]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.6
[0.1.5]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.5
[0.1.4]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.4
[0.1.3]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.3
[0.1.2]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/brainandforce/NormalForms.jl/releases/tag/v0.1.0

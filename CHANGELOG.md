# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

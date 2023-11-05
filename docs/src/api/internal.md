# Internal APIs

```@meta
CurrentModule = NormalForms
```

Most users should not need to read this section: it's primarily intended for those who would like
to contribute to the package.

Docstrings for internal functions and other objects begin with `NormalForms`, in contrast with the
public API.

## Factorizations
```@docs
NormalForms.hnf_ma!
NormalForms.snf_ma!
```

# Algorithms
```@docs
NormalForms.eye
NormalForms.gcd_kb
NormalForms.detb!
NormalForms.detb
NormalForms.find_pivot
NormalForms.is_row_zero_after
NormalForms.is_col_zero_after
NormalForms.zero_row!
NormalForms.zero_col!
NormalForms.reduce_cols_off_diagonal!
NormalForms.zero_row_and_col!
NormalForms.enforce_divisibility!
```

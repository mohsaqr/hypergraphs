# Print method for hypergraphs_bootstrap

Print method for hypergraphs_bootstrap

## Usage

``` r
# S3 method for class 'hypergraphs_bootstrap'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `hypergraphs_bootstrap` object.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  Additional arguments (ignored).

## Value

The input `x`, invisibly.

## Examples

``` r
seqs <- list(
  c("a", "b", "c", "a", "b", "c"),
  c("x", "b", "d", "x", "b", "d"),
  c("a", "b", "c", "a", "b", "c"),
  c("x", "b", "d", "x", "b", "d")
)
bs <- hg_bootstrap(seqs, n_boot = 50, max_order = 2, seed = 1)
print(bs, n = 5)
#> Memory-network bootstrap: 8 rule edges (2 higher-order) from 4 sequences
#>   50 replicates, 95% percentile CIs
#>   Higher-order rule support: min 0.54, median 0.58, max 0.62
#>  from to order count probability ci_lower ci_upper support n_boot_used
#>     a  b     1     4         1.0  1.00000  1.00000    0.94          47
#>     b  c     1     4         0.5  0.00000  0.94375    1.00          50
#>     b  d     1     4         0.5  0.05625  1.00000    1.00          50
#>     c  a     1     2         1.0  1.00000  1.00000    0.94          47
#>     d  x     1     2         1.0  1.00000  1.00000    0.96          48
#> ... 3 more rows
```

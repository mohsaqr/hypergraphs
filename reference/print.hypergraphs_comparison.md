# Print method for hypergraphs_comparison

Print method for hypergraphs_comparison

## Usage

``` r
# S3 method for class 'hypergraphs_comparison'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `hypergraphs_comparison` object.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  Additional arguments (ignored).

## Value

The input `x`, invisibly.

## Examples

``` r
set.seed(1)
first_order <- replicate(6, sample(c("a", "b", "c"), 12, replace = TRUE),
                         simplify = FALSE)
second_order <- replicate(6, rep(c("a", "b", "c", "b"), 3),
                          simplify = FALSE)
by_group <- hon(c(first_order, second_order),
                group = rep(c("first", "second"), each = 6), max_order = 2)
comparison <- hg_compare(by_group, n_perm = 99, seed = 1)
print(comparison, n = 5)
#> Memory-network comparison: first (6 sequences) vs second (6 sequences)
#>   12 pooled rule edges, 99 permutations
#>   Global weighted |diff|: 0.4814, p = 0.01 (12 of 12 rules comparable)
#>   Significant edges (BH, alpha = 0.05): 11
#>  from to order count count_first count_second probability_first
#>     a  a     1     7           7            0         0.3500000
#>     a  b     1    23           5           18         0.2500000
#>     a  c     1     8           8            0         0.4000000
#>     b  a     1    18           6           12         0.2307692
#>     b  b     1    14          14            0         0.5384615
#>  probability_second       diff p_value      p_adj significant n_perm_used
#>                 0.0  0.3500000    0.02 0.02400000        TRUE          99
#>                 1.0 -0.7500000    0.01 0.01333333        TRUE          99
#>                 0.0  0.4000000    0.01 0.01333333        TRUE          99
#>                 0.4 -0.1692308    0.07 0.07000000       FALSE          99
#>                 0.0  0.5384615    0.01 0.01333333        TRUE          99
#> ... 7 more rows
```

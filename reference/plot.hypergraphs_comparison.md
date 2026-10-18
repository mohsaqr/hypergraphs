# Plot a memory-network comparison

Difference plot of the per-edge probability differences, the largest
absolute differences first. Significance (BH-adjusted) is encoded by
both colour and point shape.

## Usage

``` r
# S3 method for class 'hypergraphs_comparison'
plot(x, top = 20L, ...)
```

## Arguments

- x:

  A `hypergraphs_comparison` object.

- top:

  Integer. Number of edges to show (by descending absolute difference).
  Default `20`.

- ...:

  Additional arguments (ignored).

## Value

A ggplot object, invisibly.

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
plot(comparison, top = 10)
```

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

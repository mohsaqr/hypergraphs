# Plot a memory-network bootstrap

Forest plot of the rule-edge probabilities with their bootstrap
percentile intervals, the most frequent edges first. Order is encoded by
both colour (Okabe-Ito) and point shape.

## Usage

``` r
# S3 method for class 'hypergraphs_bootstrap'
plot(x, top = 20L, ...)
```

## Arguments

- x:

  A `hypergraphs_bootstrap` object.

- top:

  Integer. Number of edges to show (by descending count). Default `20`.

- ...:

  Additional arguments (ignored).

## Value

A ggplot object, invisibly.

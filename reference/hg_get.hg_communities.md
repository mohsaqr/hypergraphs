# Tables of a hypergraph community ensemble

Reads, prints and plots the result of
[`hg_communities()`](https://pak.dynasite.org/hypergraphs/reference/hg_communities.md)
on a hypergraph.

## Usage

``` r
# S3 method for class 'hg_communities'
print(x, n = 10L, ...)

# S3 method for class 'hg_communities'
hg_get(
  x,
  what = c("medoid", "partitions", "runs", "sizes", "ami", "ari", "nmi", "weights"),
  ...,
  converged = NULL,
  sort_by = NULL,
  top = NULL
)

# S3 method for class 'hg_communities'
plot(x, ...)
```

## Arguments

- x:

  An `hg_communities` object.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  additional arguments passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html);
  otherwise unused.

- what:

  Table to return: `"medoid"` (default: one row per node with its
  community in the AMI-medoid run), `"partitions"`, `"runs"`, `"sizes"`,
  `"ami"`, `"ari"`, `"nmi"`, or (IRMM fits only) `"weights"`. The three
  similarity tables have one row per distinct pair of runs (`run_a`,
  `run_b`, and the similarity), without the diagonal.

- converged:

  For `what = "runs"` of an IRMM fit: `TRUE` keeps the runs whose
  weights settled, `FALSE` the runs that reached `max_iter`. `NULL`
  (default) keeps every run.

- sort_by:

  For `what = "runs"`: a numeric column of the runs table
  (`"modularity"`, `"n_communities"` or `"iterations"`) to order the
  runs by, largest first, ties broken by run number. `NULL` (default)
  keeps run order.

- top:

  `NULL` (default, every row) or the number of first rows of any table
  to return, applied after `converged` and `sort_by`.

## Value

[`hg_get()`](https://pak.dynasite.org/hypergraphs/reference/hg_get.md):
a base data.frame. Raises `hypergraphs_bad_input` for `converged` or
`sort_by` with a table other than `"runs"`, or `converged` with an
Infomap fit. [`print()`](https://rdrr.io/r/base/print.html): `x`,
invisibly. [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
cograph plot of the projection, coloured by the medoid communities.

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
  edge = rep(paste0("e", 1:4), each = 3)
)
h <- group_hypergraph(dat, "member", "edge")
if (requireNamespace("igraph", quietly = TRUE)) {
  fit <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
  hg_get(fit, what = "runs")
  hg_get(fit, what = "runs", sort_by = "n_communities", top = 1)
}
#>   run seed n_communities codelength
#> 1   1    1             2   2.052397
```

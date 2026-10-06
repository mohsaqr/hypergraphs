# Random hypergraphs

Samples a hypergraph from one of four random models, chosen with `type`.

## Usage

``` r
random_hypergraph(
  type = c("uniform", "regular", "gnp", "sbm"),
  ...,
  seed = NULL
)
```

## Arguments

- type:

  The model: `"uniform"` (default), `"regular"`, `"gnp"` or `"sbm"`.

- ...:

  The arguments of the model, by name:

  `n`

  :   Number of nodes (all types; for `"sbm"`, `NULL` uses
      `sum(block_sizes)`).

  `m`

  :   Number of hyperedges (`"uniform"`, `"regular"`, `"gnp"`; for
      `"gnp"`, `NULL` draws it from a Poisson distribution).

  `k`

  :   Hyperedge size (`"uniform"`) or node degree (`"regular"`).

  `prob`

  :   Sampling weights, length `n` for `"uniform"` and length `m` for
      `"regular"`; `NULL` is equal weights.

  `p`

  :   Incidence probability in `[0, 1]` (`"gnp"`).

  `lambda`

  :   Poisson mean of `m` (`"gnp"`).

  `P`

  :   Symmetric block-to-block pair probability matrix (`"sbm"`).

  `block_sizes`

  :   Positive integer block sizes (`"sbm"`).

  `d`

  :   Hyperedge size, recycled across sampled pairs; with
      `variable_size = TRUE`, the Poisson mean of the size minus two
      (`"sbm"`).

  `impurity`

  :   Number of added members replaced by nodes outside the blocks of
      the pair (`"sbm"`, default 0).

  `variable_size`

  :   Draw sizes as `2 + Poisson(d)` (`"sbm"`, default `FALSE`).

  `absolute_purity`

  :   Replacements come only from outside the blocks of the pair
      (`"sbm"`, default `TRUE`).

- seed:

  Optional seed; the random number state of the caller is restored on
  exit.

## Value

A `net_hg` with binary incidence and the model parameters in `$params`;
for `"sbm"`, `$blocks` records the planted block of each node.

## Details

`"gnp"` is the Erdos–Renyi-style model: `n` nodes and `m` hyperedges,
with every incidence drawn independently as Bernoulli(`p`). If `m` is
omitted it is Poisson with mean `lambda`, or `n * p` when `lambda` is
also omitted. Empty and singleton hyperedges are retained, since they
are outcomes of the model.

`"uniform"` makes every hyperedge contain exactly `k` distinct nodes,
and `"regular"` makes every node belong to exactly `k` distinct
hyperedges. Sampling weights `prob` can be unequal.

`"sbm"` samples a graph stochastic block model and extends every sampled
pair of nodes to a hyperedge with further members from the blocks of the
pair; `impurity` of the added members can then be replaced by nodes
outside those blocks.

Each type takes its own arguments, given by name; an argument that the
chosen type does not take is an error.

## Conditions

`hypergraphs_bad_input` for an argument the chosen type does not take,
or an argument given without a name.

## References

Marchette, D. J. (2021). HyperG: Hypergraphs in R. R package version
1.0.0.

## Examples

``` r
uniform <- random_hypergraph("uniform", n = 20, m = 8, k = 3, seed = 1)
hg_measures(uniform, what = "distribution", measure = "size")
#>   value n proportion ccdf
#> 1     3 8          1    1
regular <- random_hypergraph("regular", n = 20, m = 8, k = 2, seed = 1)
hg_measures(regular, what = "distribution", measure = "hyperdegree")
#>   value  n proportion ccdf
#> 1     2 20          1    1
bernoulli <- random_hypergraph("gnp", n = 20, m = 8, p = 0.2, seed = 1)
bernoulli
#> Hypergraph: 20 nodes, 8 hyperedges (sizes 1: 1, 2: 1, 3: 4, 4: 1, 5: 1)
#> Source: random gnp model
#>  hyperedge size              members
#>         h1    2              V6, V16
#>         h2    1                   V3
#>         h3    3         V5, V13, V18
#>         h4    3          V1, V7, V10
#>         h5    5 V3, V4, V8, V10, V14
#>         h6    4     V1, V9, V12, V19
#>         h7    3         V1, V14, V17
#>         h8    3         V9, V10, V13
blocks <- matrix(c(0.5, 0.05, 0.05, 0.5), 2, 2)
planted <- random_hypergraph("sbm", P = blocks, block_sizes = c(10, 10),
                             d = 3, seed = 1)
hg_get(planted, what = "nodes")
#>    node degree block
#> 1    V1      9     1
#> 2    V2      9     1
#> 3    V3     12     1
#> 4    V4      5     1
#> 5    V5      5     1
#> 6    V6      5     1
#> 7    V7      8     1
#> 8    V8      6     1
#> 9    V9      8     1
#> 10  V10      5     1
#> 11  V11      5     2
#> 12  V12      5     2
#> 13  V13      6     2
#> 14  V14     10     2
#> 15  V15     11     2
#> 16  V16     11     2
#> 17  V17      9     2
#> 18  V18      7     2
#> 19  V19     10     2
#> 20  V20     10     2
```

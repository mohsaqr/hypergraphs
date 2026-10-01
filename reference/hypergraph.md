# Build a hypergraph from any input

`hypergraph()` is the main constructor of the package. It reads the
input and builds the hypergraph the input describes, through the
constructor of that kind of data:

## Usage

``` r
hypergraph(data, ...)
```

## Arguments

- data:

  The input: a data frame, a list of sequences, a network, a topic model
  or a clustering of sequences.

- ...:

  Arguments of the constructor that the input selects.

## Value

A `net_hg` object (a `net_temporal_hypergraph` for data with a clock),
identical to the result of the constructor called directly.

## Details

- A data frame in long format with `actor` and `group` (or `from` and
  `to`) describes observed groups, and every group becomes a hyperedge
  ([`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md)).
  With `by` or `top` the groups are read as sets and the most frequent
  sets become hyperedges. With `time`, `start` or `end` the hyperedges
  carry a clock
  ([`temporal_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/temporal_hypergraph.md)).
  With `window`, `step` or `action` the data are sequences, and every
  window of consecutive actions becomes a hyperedge
  ([`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)).

- A list of sequences is read the same way
  ([`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)).

- A network, as a weight matrix, a sparse matrix, a `netobject` or a
  `cograph_network`, has its cliques promoted to hyperedges
  ([`network_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/network_hypergraph.md));
  with `window`, a model object built from sequences is read as
  sequences instead.

- A topic model fitted by
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md),
  or a clustering of sequences, becomes the hypergraph of its frequent
  sets
  ([`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md)).

Every argument in `...` is passed to that constructor, whose
documentation describes it, and the result is the constructor's own.

## Conditions

`hypergraphs_bad_input` for an input of a class no constructor reads,
naming the class, and the conditions of the constructor called.

## References

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
[doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004)

## See also

[`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md),
[`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md),
[`temporal_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/temporal_hypergraph.md),
[`network_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/network_hypergraph.md),
[`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md)
for a corpus.

## Examples

``` r
meetings <- data.frame(
  person = c("Alice", "Bob", "Carol", "Alice", "Bob", "Dave", "Eve"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m3", "m3"))
hypergraph(meetings, actor = "person", group = "meeting")
#> Hypergraph: 5 nodes, 3 hyperedges (sizes 2: 2, 3: 1)
#> Source: group membership (actor = person, group = meeting)
#>  hyperedge size           members weight
#>         m1    3 Alice, Bob, Carol     NA
#>         m2    2        Alice, Bob     NA
#>         m3    2         Dave, Eve     NA

sessions <- list(c("a", "b", "c", "a"), c("b", "c", "d"))
hypergraph(sessions, window = 2)
#> Hypergraph: 4 nodes, 4 hyperedges (sizes 2: 4)
#> Source: windowed sequences, window = 2, step = 1 (5 windows from 2 sequences)
#>  hyperedge size members weight
#>         h1    2    a, b      1
#>         h2    2    a, c      1
#>         h3    2    b, c      2
#>         h4    2    c, d      1

weights <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3,
                  dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
hypergraph(weights)
#> Hypergraph: 3 nodes, 4 hyperedges (sizes 2: 3, 3: 1)
#> Source: network cliques (p = 1.00, include_pairwise = TRUE, max_size = 3)
#>  hyperedge size members weight
#>         h1    2    y, z     NA
#>         h2    2    x, y     NA
#>         h3    2    x, z     NA
#>         h4    3 x, y, z     NA
```

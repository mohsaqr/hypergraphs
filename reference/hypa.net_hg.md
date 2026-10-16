# Hypergeometric anomaly detection for co-occurring node pairs

The hypergraph method of
[`hypa()`](https://pak.dynasite.org/hypergraphs/reference/hypa.md).
Scores pairs of nodes against a hypergeometric null that keeps every
node's co-occurrence strength and otherwise pairs the nodes at random. A
pair that shares far more hyperedges than the null predicts is
over-represented; a pair that shares far fewer, possibly none, is
under-represented. This is the undirected counterpart of
[`hypa.default()`](https://pak.dynasite.org/hypergraphs/reference/hypa.default.md),
which applies the same null to the paths of a higher-order network, and
it reports the same columns so the two read alike.

## Usage

``` r
# S3 method for class 'net_hg'
hypa(x, min_count = 5L, alpha = 0.05, p_adjust = "none", top = NULL, ...)
```

## Arguments

- x:

  A `net_hg`, or a snapshot of a temporal one.

- min_count:

  Single whole number \\\ge 1\\ (default `5L`). A pair is
  over-represented only if it co-occurs at least `min_count` times and
  under-represented only if it is expected to co-occur at least
  `min_count` times. Pairs eligible for neither are not scored.

- alpha:

  Significance threshold for the `anomaly` column, a single number in
  `(0, 0.5)` (default `0.05`), applied to the adjusted p-values.

- p_adjust:

  Multiplicity correction, any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method or
  `"none"` (default `"none"`, as in LaRock et al. 2020). `"BH"` controls
  the false discovery rate within each direction.

- top:

  Optional number of rows to keep, the most significantly
  over-represented first. `NULL` (default) keeps every scored pair.

- ...:

  Must be empty: an argument that only the sequence method takes
  (`order`, `action`, ...) raises `hypergraphs_bad_input`.

## Value

A base `data.frame`, one row per scored pair, most significantly
over-represented first (by `p_adjusted_over`, ties broken by `ratio`,
then by name), with columns `from`, `to` (the pair, in the node order of
the hypergraph), `observed` (co-occurrence count, `0` for a pair that
never co-occurs), `expected` (the fitted \\E\_{ij}\\), `ratio` (observed
/ expected), `p_under`, `p_over`, `p_adjusted_under`, `p_adjusted_over`
and `anomaly` (`"over"`, `"under"` or `"normal"`). The attribute
`"hypa_fit"` records the fit of the null: `converged`, `iterations`,
`max_error` (largest absolute strength error) and `n_cooccurrences`
(\\m\\).

Errors of class `hypergraphs_bad_input` are raised for an invalid
argument, and of class `hypergraphs_empty_result` when no pair co-occurs
or no pair reaches `min_count` in either sense. A warning of class
`hypergraphs_no_converge` is raised when the null cannot be fitted to
tolerance.

## Details

The test is analytic, so it needs no resampling. That is what makes it
usable where
[`hg_null_test()`](https://pak.dynasite.org/hypergraphs/reference/hg_null_test.md)
is not:
[`hg_null_test()`](https://pak.dynasite.org/hypergraphs/reference/hg_null_test.md)
answers whether the hypergraph as a whole carries more repetition than
chance, by Monte Carlo, and returns one row; this returns one row per
pair.

**Observed counts.** For nodes \\i \< j\\, \\c\_{ij}\\ is the number of
hyperedges that contain both (hyperedge weights are ignored; a repeated
hyperedge counts each time). The total is \\m = \sum\_{i\<j} c\_{ij}\\
and the strength of node \\i\\ is \\s_i = \sum\_{j \ne i} c\_{ij}\\, so
\\\sum_i s_i = 2m\\.

**Null model.** Every unordered pair of distinct nodes with positive
strength is a possible cell, whether it co-occurs or not. The expected
count is \\E\_{ij} = x_i x_j\\, the product form of Chung and Lu (2002),
with the fitness \\x\\ fitted so that \\\sum\_{j \ne i} E\_{ij} = s_i\\
for every node: the null preserves each node's observed strength, as the
propensity fit of LaRock et al. (2020) preserves the in- and
out-strengths of a De Bruijn graph. The fit is the damped fixed point
\\x_i \leftarrow \sqrt{x_i s_i / (X - x_i)}\\, \\X = \sum_k x_k\\,
iterated until the largest strength error falls below \\10^{-10} m\\ or
10,000 sweeps have run. Each sweep costs \\O(n)\\; no \\n \times n\\
matrix is formed. With three or more nodes the fit has a solution with
positive fitness exactly when no node takes part in every co-occurring
pair. In a star-like hypergraph, where one node is in every co-occurring
pair, the fit cannot converge, and the result is returned with a warning
of class `hypergraphs_no_converge`.

**Test.** The propensity of a pair is \\\Xi\_{ij} = m E\_{ij}\\, which
sums to \\m^2\\. The observed count \\f = c\_{ij}\\ is referred to a
\\\mathrm{Hypergeometric}(N, K, n)\\ law with \\N = m^2\\, \\K =
\mathrm{round}(\Xi\_{ij})\\ and \\n = m\\ draws, whose mean is
\\E\_{ij}\\. `p_under` is \\P(X \le f)\\ and `p_over` is the strict
upper tail \\P(X \> f) = 1 - \\ `p_under`, so a pair is over-represented
when its HYPA score \\P(X \le f)\\ exceeds \\1 - \alpha\\, the rule of
LaRock et al. (2020). With small counts the strict upper tail rejects
more often than \\\alpha\\ under the null; `min_count` limits the
problem.

**Scored pairs and classification.** A pair is scored when it could be
classified: when it co-occurs at least `min_count` times (eligible for
over-representation) or is expected to co-occur at least `min_count`
times (eligible for under-representation, including pairs that never
co-occur). There are at most \\m /\\ `min_count` pairs of the second
kind, because the expected counts sum to \\m\\. A scored pair is
`"under"` when `expected >= min_count` and `p_adjusted_under < alpha`,
otherwise `"over"` when `observed >= min_count` and
`p_adjusted_over < alpha`, otherwise `"normal"`. `min_count` selects the
hypotheses by count, never by p-value, and the null is fitted on every
pair whatever the floor, so raising the floor changes which pairs are
tested but never the expected counts.

**Multiple testing.** LaRock et al. (2020) classify against `alpha`
without correction, and so does the default `p_adjust = "none"`. Any
[stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method can be
applied instead. The under- and over-representation p-values are then
adjusted separately, and each family is every possible pair
(\\n(n-1)/2\\ for the \\n\\ nodes with positive strength), as in the
sequence method, which adjusts over every possible path. Only the scored
pairs are computed; the others enter
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) through its
`n` argument, which treats them as larger p-values. Adjusting within the
scored pairs alone would under-correct, because a pair enters the scored
set partly through its observed count.

Time is handled by scoring a cumulative snapshot rather than by slicing
into periods. Disjoint periods split a pair's evidence across cells
while multiplying the number of tests, so both terms move the wrong way;
a sequence of `hg_snapshot(hg, at = ..., mode = "cumulative")`
accumulates evidence instead and shows when a pair becomes anomalous.

## References

LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies
in time series data on networks. *Proceedings of the 2020 SIAM
International Conference on Data Mining*, 460-468.
[doi:10.1137/1.9781611976236.52](https://doi.org/10.1137/1.9781611976236.52)

Chung, F., and Lu, L. (2002). Connected components in random graphs with
given expected degree sequences. *Annals of Combinatorics*, 6(2),
125-145. [doi:10.1007/PL00012580](https://doi.org/10.1007/PL00012580)

## See also

[`hypa.default()`](https://pak.dynasite.org/hypergraphs/reference/hypa.default.md)
for the same null on the paths of sequences,
[`hg_null_test()`](https://pak.dynasite.org/hypergraphs/reference/hg_null_test.md)
for the Monte Carlo whole-hypergraph tests.

## Examples

``` r
memberships <- data.frame(
  actor = c("a","b","c", "a","b","d", "a","b","e", "c","d","e", "a","b","f"),
  group = rep(c("g1","g2","g3","g4","g5"), each = 3)
)
hg <- group_hypergraph(memberships, node = "actor", hyperedge = "group")
hypa(hg, min_count = 2L)
#>   from to observed expected    ratio   p_under    p_over p_adjusted_under
#> 1    a  b        4 3.354836 1.192309 0.7801609 0.2198391        0.7801609
#>   p_adjusted_over anomaly
#> 1       0.2198391  normal
```

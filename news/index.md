# Changelog

## hypergraphs 0.6.1

- **The package is renamed from hypernets to hypergraphs.** Every verb
  keeps its name; condition classes and result classes change their
  prefix from `hypernets_` to `hypergraphs_` (`hypergraphs_bad_input`,
  `hypergraphs_result`, …), so code that catches a condition by class
  changes with them.

- **[`hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/hypergraph.md)
  is the main constructor.** It reads any input and calls the
  constructor of that kind of data, returning its result unchanged: a
  data frame of groups goes to
  [`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md)
  (with `by` or `top`, the frequent sets), with a clock to
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/temporal_hypergraph.md),
  with `window`, `step` or `action` to
  [`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md);
  a list of sequences to
  [`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md);
  a network (matrix, sparse matrix, `netobject`, `cograph_network`) to
  [`network_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/network_hypergraph.md);
  a topic model or a clustering of sequences to
  [`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md).
  The specific constructors remain.

- **`group_hypergraph(by =)` counts frequent sets within groups.** With
  `by`, each value of `group` (a trial, a session, a basket) is one set
  of `actor` values, and the `top` most frequent sets within each value
  of `by` (an outcome, a cluster) become hyperedges, read with
  `hg_get(x, what = "sets")` and drawn per value with
  `plot(x, group = )`, as for clustered sequences. The topic-combination
  route names its grouping column `by` as well. Giving `top` without
  `by` counts the sets over all groups. New dataset `debug_events`:
  4,000 sessions with a coding assistant, one row per event, with runs
  (a stretch of a session between two verdicts) and a quick or slow
  group; real event sequences of another domain with every event renamed
  (`data-raw/debug_events.R`).

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  hypergraph colours `color_by = "size"` discretely, one Okabe-Ito
  colour per size; it was a continuous scale.

- `network_hypergraph(type = "vr")` raises `hypergraphs_bad_input` and
  points to `simplicial(type = "vr")`; it raised a plain error that said
  the package builds no Vietoris-Rips complex.

- **Communities of memory networks can be drawn over the transition
  network.** [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  a
  [`hg_communities()`](https://mohsaqr.github.io/hypergraphs/reference/hg_communities.md)
  result can draw the network as a transition network analysis plot (TNA
  styling through
  [`cograph::overlay_communities()`](https://sonsoles.me/cograph/reference/overlay_communities.html)):
  arrows carry transition probabilities, node area follows flow, and
  every community is a blob around its nodes. New `type = "network"`
  draws the network of states, in which a shared state lies in two
  blobs; `type = "states"` now draws the memory nodes the same way, laid
  out by the walk’s transitions between them, in place of cograph’s
  group-ring layout. The default pebble view of the community hypergraph
  is unchanged, and the labels say “state” and “memory node” as the
  tables do.

- [`hg_network()`](https://mohsaqr.github.io/hypergraphs/reference/hg_network.md)
  replaces
  [`hg_relations()`](https://mohsaqr.github.io/hypergraphs/reference/hg_relations.md)
  (kept as a deprecated alias that warns with `hypergraphs_deprecated`
  and returns the identical result). Besides a partition (`clusters =`),
  it takes a topic model (`topics =`). Without `threshold` two topics
  are linked by the correlation of their shares (stm’s simple
  `topicCorr()`, equal on the local oracle); with `threshold` the weight
  is the number of documents in which both topics reach that share
  (Abuhay et al. 2017; Cassi et al. 2017), normalisable by the same
  similarity measures.

- `group_hypergraph(topic_model, threshold =)` builds the hypergraph of
  topic combinations: each document is the set of its topics at the
  threshold, and the most frequent sets become hyperedges, counted as
  for clustered sequences (`hg_get(x, what = "sets")`). `by =` counts
  them within a column of the documents table (a year, an author);
  `min_size` and `top = Inf` keep the multi-topic sets and all of them.
  Hypergraphs of clustered sequences now print their source.

- `hg_sequences(topics =)` builds sequences from each document’s main
  topic, grouped by any column of the documents table.

- [`hg_markov_stability()`](https://mohsaqr.github.io/hypergraphs/reference/hg_markov_stability.md)
  reads only the possible transitions of a matrix when it checks
  irreducibility, so an unnormalised matrix warns once with
  `normalize = TRUE` (it warned twice) and is refused without a warning
  with `normalize = FALSE` (it warned that it was normalising).

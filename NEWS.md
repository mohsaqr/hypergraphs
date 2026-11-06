# hypergraphs (development version)

## Breaking changes

* Column names give no role, as in `Nestimate::build_network()` since Nestimate
  1.0.7, so the same arguments build the same sequences in both packages. A
  column called `action`, `time`, `session` or `session_id` was used when its
  argument was `NULL`. Now a column called `action` that is not declared raises
  `hypergraphs_long_format`, and an undeclared `time`, `session` or `session_id`
  column beside `action` is named in a message of class
  `hypergraphs_column_role` and not used: events follow row order and sessions
  are not split. Calls that relied on those names need the arguments added.
  `session = FALSE` now means the same as `NULL`. Requires Nestimate >= 1.0.7.
* `temporal_hypergraph()`, `group_hypergraph()` and `hypergraph()` on a table no
  longer detect columns from an alias table (`source`/`target`, `Timestamp`,
  `onset`/`terminus`, `actor`/`group`, a `weight` column, ...). Name the columns:
  `from` and `to`, or `node` and `hyperedge`, and `time`, or `start` and `end`.
  A relational table without them raises `hypergraphs_bad_input`, also when it
  is all categorical and would otherwise read as sequences. A `weight` column
  counts only when named. A `nodes` table gives entry times in its `start`
  column; `time` and `date` are no longer read.
* `hg_centrality()` now defaults to `type = "clique"`. The `"Z"` and `"H"`
  centralities are defined by Benson (2019) for uniform hypergraphs only:
  with mixed hyperedge sizes the equation mixes polynomial degrees, so its
  solution depends on the normalisation chosen inside the iteration.
  Requesting them on a non-uniform hypergraph raises
  `hypergraphs_nonuniform_hypergraph`; `hg_keywords(centrality = "Z" | "H")`
  does the same for non-uniform word hypergraphs. The main vignette now
  reports the clique centrality and PageRank for the trial hypergraph,
  whose trials hold two to six events.
* The `"Z"` and `"H"` iterations normalise the tensor update before the
  stabilising shift and require a small eigen-equation residual to count as
  converged. On large uniform hyperedges the shift swamped the update and
  returned the uniform vector without a warning (a 20-member, 3-petal
  sunflower gave a hub/petal ratio of 1 instead of sqrt(3)).
* The memory and topology verbs follow the Nestimate corrections of
  2026-10-27 once that version is installed: path-based degrees of freedom
  in `mogen()` and `markov_order()`, the Scholtes order rule for
  `mogen(criterion = "lrt")`, Atkin's Q-vector in `hg_qanalysis()`, and
  persistent homology built one dimension above `max_dim`.

* A function-by-function audit (2026-11-01) tightened every contract.
  Arguments a verb could not use were silently ignored before and are now
  refused with `hypergraphs_bad_input`: `simplicial()` controls that its
  `type` does not read, `group_hypergraph()` options on clustering or
  topic input, `hg_bootstrap()` estimation arguments on a grouped model,
  `hg_get()` filters that the chosen table does not read (`order_min`,
  `sort_by`, `topic`, `dimension`, `normalized`, ...), unused plot
  arguments, and a `type` named twice. Every `seed` is validated the same
  way (one whole number) and restored through one helper. A long table
  found by its `action` column but holding an unnamed actor column now
  raises `hypergraphs_long_format` instead of joining every actor into one
  sequence, matching the rule that a long table names its roles. Text that
  is not valid UTF-8 is refused by name instead of crashing inside
  `tolower()`. `hg_agreement()` returns `NA` with
  `hypergraphs_undefined_statistic` when the partitions share one node.
  `clean_text()` cuts only a real copyright notice (a mark followed by a
  year, a copyright sign or "the author(s)"), the last one, so "(c)" list
  items and sentences about copyright law survive; 7 of the 1000
  `covid_sample` abstracts change by a stray sign.

## Bug fixes

* Wrong answers found by the audit: `window_hypergraph()` split a state
  label containing the internal separator (such labels are now refused);
  an edge list row with one missing end became a singleton hyperedge in
  `group_hypergraph()` and `temporal_hypergraph()` (the row is dropped with
  `hypergraphs_dropped_rows`); `hg_keywords()` on a kNN text hypergraph
  ranked document ids as words; `hg_centrality(type = c("clique",
  "clique"))` returned a duplicated column; the sparse pair reader of
  `pairwise_network()` lost every pair of a lower-triangle symmetric
  matrix; `hg_pagerank()` refused a vertex in no hyperedge where
  `hg_centrality(type = "pagerank")` ranked it (dangling mass now follows
  the teleport vector, so the two agree); `hg_topic_search(parallel =
  TRUE)` lost the `hypergraphs_no_converge` warning raised in a worker;
  `clean_text(min_chars =)` and its stop-word step cut "p53" to "53";
  `hg_markov_stability()` refused a grouped model it documents; `hon()`,
  `mogen()` and `hypa()` read a `netobject_group` or any classed list as
  sequences and fitted a network over garbage states; `hg_get(mogen_fit,
  what = "transitions")` and `summary()` crashed when order 0 was selected;
  `summary()` on an `hg_motifs()` table crashed in `print.default`;
  `plot(what = "sizes")` of a community comparison drew nothing for a fit
  whose communities share one size; `plot()` with `top = 0` crashed on
  bootstrap and comparison results; `hg_homology(type = "rips")` failed
  although documented.
* New where the audit found a gap: `plot()` for a grouped bootstrap
  (`group =`), `hg_get()` for a grouped Markov stability result,
  `hypergraph()` on a `tna` model, `print()` for a motif test table,
  `hg_get(top =)` on temporal hypergraphs and result summaries, and a
  classed `hypergraphs_empty_result` from `window_hypergraph()` and
  `hg_line_graph()` on empty input.

* Memory constructors preserve literal numeric states from Nestimate models.
  Only tna models and explicitly encoded converted data decode positional
  state codes. Equivalent raw and model inputs now give the same HON,
  MOGen and HYPA results.
* Model input is normalized locally before calling Nestimate, so this fix
  also works with the minimum supported Nestimate 0.9.26.
* Grouped mixture/clustering input preserves numeric state values.
  Component matrix labels alone no longer imply encoded data; explicit
  `labels` or `alphabet` attributes identify positional codes.

# hypergraphs 0.7.10

First release on CRAN.

hypergraphs analyzes higher-order networks, in which a relation binds more
than two nodes or depends on more than the current node. It covers four
families of structures.

* Memory networks from categorical sequences: higher-order networks with
  rule extraction (`hon()`), their embedding (`honem()`), the multi-order
  generative model (`mogen()`), path anomaly detection (`hypa()`), a
  permutation test of the Markov order (`markov_order()`), centralities and
  communities of a higher-order network, and bootstrap and permutation
  inference for its rules (`hg_bootstrap()`, `hg_compare()`).
* Simplicial complexes: clique, Vietoris-Rips, window and pathway complexes
  (`simplicial()`), Betti numbers, persistent homology, persistence
  landscapes, bottleneck and Wasserstein distances, and q-analysis.
* Hypergraphs: construction from group memberships, sliding windows,
  network cliques, interval records and embedding neighbors
  (`hypergraph()` and the `*_hypergraph()` constructors); structural
  measures, centralities, Laplacians, spectral and symmetric-NMF
  clustering, label spreading, modularity, mixed-membership communities,
  motifs, null models and temporal snapshots.
* Text hypergraphs: a corpus as a document-word hypergraph
  (`text_hypergraph()`, `clean_text()`), analyzed by clustering, keyword
  extraction, topic models (`hg_topics()`) and classification from a few
  labels. A clustered corpus becomes a set of sequences for the memory
  family with `hg_sequences()`.

Every result is read with `hg_get()`, which returns a data frame, and every
function that takes sequences accepts long, wide, list or model input.

# hypergraphs (development version)

## Breaking changes

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

## Bug fixes

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

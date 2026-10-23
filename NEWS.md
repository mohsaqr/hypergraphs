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

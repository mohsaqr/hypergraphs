# hypergraphs 0.7.2

* The package has no compiled code again: `hg_topics()` computes its
  expected counts in R, with the same values as 0.7.1.

# hypergraphs 0.7.1

Corrections from a function-by-function audit. Results change only where
they were wrong; valid input that was handled correctly gives identical
output.

## Speed

* `hg_communities()` is faster with identical results. IRMM builds each
  reweighting pass's clique reduction as a sparse matrix (about 3 times
  faster on `icsid_tribunals`), and the AMI between runs computes the
  expected mutual information in one vectorised pass. New `parallel` and
  `n_cores` arguments run the independent seeded runs with
  `parallel::mclapply()` (not on Windows); the fit is identical to the
  serial one. On `icsid_tribunals` (441 arbitrators, 742 cases), 8 cores:
  Infomap 34 s to 11 s, IRMM 27 s to 3 s.
* `hg_mmsbm()` gains `parallel` and `n_cores`: every start's initial
  values are drawn first, in order, from one stream, then the EM fits run
  with `parallel::mclapply()`. Fits and the caller's random stream are
  identical to the serial run. `hg_topics()` runs its parallel starts
  through the same helper, so a failed worker raises
  `hypergraphs_parallel_failed` instead of returning an error object as a
  start.
* `hg_topics()` and `hg_topic_search()` no longer exhaust memory on a large
  corpus. The expected counts are gathered in blocks of about 80 MB
  (identical values); a 6,630-document, 91,755-word corpus at `k = 52`
  allocated about 8 GB per step before. Topic quality keeps the document x
  word presence sparse and densifies only the top words (it was a dense
  4.9 GB matrix on that corpus), and matching topics across starts is
  computed for all pairs at once (0.53 s to 0.005 s per comparison at
  `k = 52`, identical values).
* `hg_cluster(algorithm = "symnmf")` and `hg_embed(method = "symnmf")` run
  on a sparse hypergraph instead of refusing it with
  `hypergraphs_dense_required`: the dense similarity is built from the
  sparse Laplacian (up to 15,000 nodes; beyond that
  `hypergraphs_sparse_too_large`). On the corpus above, SymNMF graded
  memberships for `k = 16` take about a minute per start.
* `hg_cluster(algorithm = "symnmf")` gains `parallel` and `n_cores`: every
  start's initial factor is drawn first, so parallel starts give the serial
  result.
* k-means in `hg_cluster()` (spectral) and `hg_cocluster()` finishes a
  Hartigan-Wong start that stopped early ("Quick-TRANSfer stage steps
  exceeded" or the iteration limit) from the centres it reached, instead
  of comparing it unfinished and warning; when every start finishes the
  result is identical to `stats::kmeans()`.
* `text_hypergraph()` keeps a data column named `doc`, `n_tokens` or
  `n_types` as metadata `input_doc`, `input_n_tokens`, `input_n_types`,
  with a `hypergraphs_renamed_column` warning, instead of refusing the
  table.
* Parallel runs on macOS: a forked worker the system kills (Apple's
  Accelerate BLAS is not fork-safe) is rerun serially with a
  `hypergraphs_parallel_fallback` warning, and the result is unchanged.
* `hg_mmsbm()` explains a collapsed membership: the rate of a hyperedge sums
  over its node pairs, so in hyperedges of three or more nodes the other
  members can explain it and a node's activity goes to zero. The warning
  and the documentation say so, and point to `hg_topics()` for document
  mixtures, instead of suggesting more starts.

## Package

* The package ships one vignette, `vignette("hypergraphs")`. The text
  hypergraph guide, the COVID-19 topic walkthroughs, the text constructions
  and the document classification guide are articles on the package
  website, which keeps the source package small and its check
  short; the saved classification results moved with their article. `hg_hypergat()`,
  `hg_neural()` and the classification reader have runnable examples on
  the bundled `forum_posts`.

## Construction and temporal hypergraphs

* Hyperedges are identified by their members. Two actor/session pairs
  whose pasted labels coincide (`a.b` + `c`, `a` + `b.c`), a node named
  `a + b` beside the pair `{a, b}`, and two sequences sharing a list name
  were merged into one hyperedge. Sequence lists now need unique, non-empty
  names, and missing actor or session identifiers raise
  `hypergraphs_bad_input`.
* Metadata columns named like a structural column (`edge`, `node`,
  `start`, `weight`, ...) no longer replace the selected structure; they
  are kept under a suffixed name. Edge attributes with missing values take
  the first non-missing value of the edge whatever the row order.
* Membership weights must be finite, non-negative numbers; negative,
  infinite, factor and character weights raise `hypergraphs_bad_input`.
* `network_hypergraph()` reorders a matrix whose column names are a
  permutation of its row names, and refuses differing, duplicated or
  missing labels. Column labels were previously overwritten, which moved
  edges.
* `temporal_hypergraph()` refuses clocks that fail to parse, infinite
  clocks and per-membership intervals that end before they start. ISO 8601
  offsets (`Z`, `+0200`, `+02:00`) and fractional seconds are honoured,
  and trailing text is an error (`hypergraphs_unparsed_time`). An all-missing
  `end` column means every membership is open. Open and closed memberships
  of one hyperedge are filtered by their own spells, so snapshots no longer
  depend on row order.
* `hg_growth()` and `hg_snapshot()` share one window boundary, the
  temporal `summary()` counts distinct members, and `hg_edges()` and the
  readers of an empty snapshot return typed zero-row tables.
* The canonical `node` and `hyperedge` columns are detected without being
  named. `knn_hypergraph()` requires a whole `k` and finite embeddings.
* `dual_hypergraph()` is the transpose of the full incidence, isolated
  vertices and empty hyperedges included, and a sparse dual is a complete
  `net_hg`. `hg_subset()` keeps planted SBM blocks aligned with the kept
  nodes. `summary()` and `hg_laplacian()` accept sparse hypergraphs.
* `plot(pieces = "row")` places isolated nodes.

## Measures, random models and null tests

* Sparse `hg_measures()` equals the dense result (pairwise participation,
  isolate edge sizes, uniform density). `hg_null_test(statistic =
  "density")` uses the same density as `hg_measures()`.
* `random_hypergraph(type = "sbm")` and the configuration null drew from
  `1:x` when a candidate pool held one node, producing wrong edge sizes and
  invented memberships. Draws from larger pools are unchanged.
* `hg_null_test()` requires a whole `n`, and a statistic undefined on the
  input (`avg_jaccard` with fewer than two hyperedges) is `NA` with a
  `hypergraphs_undefined_statistic` warning.
* Clique eigenvector centrality is the Perron vector of each connected
  component; power iteration oscillated on bipartite clique expansions. Z-
  and H-eigenvector iterations that reach `max_iter` warn with
  `hypergraphs_no_converge`. Subhypergraph centrality is computed per node
  in log space, so isolated nodes are 0 rather than `-Inf`. A hypergraph of
  singleton hyperedges has zero clique, Z and H centrality.
* Empty hyperedges contribute nothing to Laplacians, random walks and
  PageRank (they produced NaN). `hg_pagerank()` refuses non-finite or
  conflicting `personalized` weights and, at `damping = 1`, disconnected
  input (`hypergraphs_hypergraph_disconnected`).
* Assortativity is computed on centred scores, so near-regular hypergraphs
  of high degree are no longer reported as undefined.
* `hg_compare_communities()` reports `NA` agreement for fits sharing fewer
  than two nodes and scores each medoid on the projection its fit saved.
  Directed citation fits have an `NA` quality row. `edge_source` is
  deprecated there. Duplicated node or source assignments raise
  `hypergraphs_bad_input` in every partition reader.
* `hg_mmsbm()` does not declare convergence while memberships collapse.
* The cluster plot supports more than nine clusters (Okabe-Ito colours
  recycle with distinct shapes), and the agreement heatmap shows negative
  AMI and ARI on a diverging scale.

## Memory networks and simplicial complexes

* Missing actions split a trajectory into contiguous runs. No transition
  crosses a gap and no state `NA` is created; `hon()`, `mogen()`,
  `markov_order()` and `hypa()` accept sequences with gaps, and
  `hg_bootstrap()` resamples whole trajectories. A real state spelled
  `"NA"` is kept.
* State labels containing ` -> ` or the internal separators raise
  `hypergraphs_bad_input`, since higher-order node names are built from
  them.
* `hg_compare()` averages over the rules whose context both groups
  observe and reports their number; with none the statistic and p-value
  are `NA` (`hypergraphs_undefined_statistic`). The global statistic of
  partially overlapping groups changes accordingly.
* Window persistence (`hg_homology()` on `simplicial(type = "window")`)
  keeps essential classes essential, and every Betti row comes from the
  same Z/2 intervals. `hg_wasserstein()` treats a finite death at 0 as
  finite outside clique mode and computes large orders without overflow.
* `simplicial(validate = TRUE)` refuses an automatic shuffle count above
  100,000 and asks for an explicit `n_null`.
* `hg_bootstrap()` with no rule above `min_freq` raises
  `hypergraphs_empty_result`. HONEM variance and Infomap savings are 0 on
  zero spectra and zero code lengths.
* `?memory` states that `order` is the number of conditioning states;
  `?hg_betti` states that Betti numbers are computed over the rationals and
  differ from `hg_homology()` on torsion.

## Text hypergraphs and neural classifiers

* `clean_text()` keeps alphanumeric tokens (`covid19`, `p53`, `covid-19`)
  intact and replaces invalid numeric entities with U+FFFD.
* `hg_keywords()` and `hg_topic_quality()` handle clusters without
  eligible words and one-word vocabularies; external keyword scores count
  distinct documents.
* `text_hypergraph()` refuses metadata columns named `doc`, `n_tokens` or
  `n_types`.
* `hg_cocluster()` checks `k` against the singular vectors available.
  `topic_network()` keeps isolated topics and treats constant topics as
  isolates (`hypergraphs_constant_topics`).
* A topic model is checked against the documents, words and counts it was
  fitted on (`net_hg_topics` gains `$corpus`). KL-NMF ignores explicit
  sparse zeros. `hg_membership()` sums to one at coincident centres.
  `hg_topic_sizes()` refuses negative weights.
* Label inputs of every classifier and topic reader refuse conflicting
  duplicates; an `NA` label marks a node as unlabelled. Embedding and
  feature row names must be unique.
* `hg_hypergat(min_count =)` counts token occurrences, as documented,
  rather than sentences.
* `hg_neural()` restores an absent random seed, HNHN operators are
  normalised in log space, and the confusion table keeps a class named
  `(unscored)`.

## Arguments

* Counts (`n`, `k`, `nstart`, `max_iter`, `n_boot`, `n_perm`, `top`,
  `dimension`, ...) must be whole numbers in range; fractions, `-Inf` and
  values beyond the integer range raise `hypergraphs_bad_input` instead of
  being truncated.

# hypergraphs 0.6.10

* The result classes of the memory family drop the `hon` prefix left over
  from the package's former name: `net_hon_boot` is now
  `hypergraphs_bootstrap`, `net_hon_boot_group` `hypergraphs_bootstrap_group`,
  `net_hon_compare` `hypergraphs_comparison`, `net_hon_communities`
  `hypergraphs_memory_communities` and `net_hon_group`
  `hypergraphs_memory_group`. Code that tests these classes with
  `inherits()` needs the new names. `hon()` and `honem()` keep their names
  (the published BuildHON and HONEM methods), and `net_hon` remains the
  Nestimate class that `hon()` returns. Printed headers read
  "Memory-network bootstrap" and "Memory-network comparison".

# hypergraphs 0.6.9

* The Argentina example article gains a storyline of the busiest
  arbitrators and restores the centrality, per-community modularity and
  growth tables of its first version.
* The repository now holds only package files; development notes and
  local tooling are kept out of it.

# hypergraphs 0.6.8

* Test suite only: the regression test of default hull plots, whose baseline
  was recorded on macOS, is skipped on Linux and Windows, where the
  force-directed layout settles about 1e-4 away.

# hypergraphs 0.6.7

* `hg_communities()` on a memory network numbers communities of equal flow
  the same way on every platform: flows equal to 12 significant digits are
  ties, broken by the first node.

# hypergraphs 0.6.6

* `window_hypergraph()` gains `collapse`: `FALSE` keeps every window as its
  own hyperedge, in order, named by its positions (`"1-3"`), with
  `sequence`, `start` and `end` in the edge metadata.
* `plot()` on a hypergraph gains `type = "storyline"` for hyperedges with an
  order (the stored order or `sort_by`), such as the windows of one
  sequence: codes become lines that the windows gather.
* `plot()` of a temporal hypergraph gains `type = "storyline"` (Tanahashi
  and Ma 2012): the `top` nodes with the most hyperedges (default 8) are
  lines from their first to their last hyperedge, hyperedges are columns in
  order of their start, and each column gathers the lines of its members.
  Lines are ordered by the barycentre rule (Sugiyama et al. 1981) over
  repeated sweeps, keeping the ordering with the fewest crossings. `start`
  and `end` limit the period; `edge_labels` and `point_size` style it.
  `spacing = "strength"` brings neighbouring lines closer the more
  hyperedges their nodes share (a full row for none, 0.35 for the strongest
  tie in the plot, linear in between) without changing their order.
  `width_by = "degree"` widens each line with its node's number of
  hyperedges in the period, with a width legend. Each
  line has an Okabe-Ito colour and a point shape, distinct for up to 72
  lines, named in the legend.
  `type` is now an explicit argument of `plot.net_temporal_hypergraph()`,
  and an argument of one view given to another raises
  `hypergraphs_bad_input`.

# hypergraphs 0.6.5

* Requires Nestimate (>= 0.9.1), whose `prepare()` takes the `timezone`
  argument the sequence input passes; with 0.8.5 the vignettes failed.

# hypergraphs 0.6.4

* `plot()` on a hypergraph gains `type = "incidence"`, the incidence matrix
  after UpSet (Lex et al. 2014): one row per node ordered by hyperdegree, one
  column per hyperedge, a bar joining each hyperedge's members, hyperdegree
  bars beside the rows and size bars above the columns when sizes differ.
  It stays legible where hulls overlap. The new `sort_by` orders the columns
  (`"size"`, an edge-metadata column such as a date, or one value per
  hyperedge; numbers largest first). `color_by`, `labels`, `edge_labels`,
  `label_size`, `edge_label_size` and `legend_title` apply; a hull-only
  argument raises `hypergraphs_bad_input`. Counts sort largest first;
  dates, text and clock columns (`start`, `end`, ...) in time order.
* `plot()` on a hypergraph sizes the node labels to the number of labelled
  nodes when `label_size` is not given: 4.2 mm up to 12 labels, then
  shrinking with the square root of the count to 2.2 mm from 44 labels on
  (was a fixed 4.2 mm). An explicit `label_size` is used as given; an
  invalid one raises `hypergraphs_bad_input`.
* `plot()` of a temporal hypergraph plots its snapshot as a hypergraph with
  `plot.net_hg()` (hulls, or `type = "incidence"`), so the hyperedges are
  kept; it used to plot a pairwise projection through `cograph::splot()`.
  `method` is deprecated (`hypergraphs_deprecated`);
  `plot(pairwise_network(hg_snapshot(x, at)))` plots the projection.
* `hg_snapshot()` keeps the data's names for nodes and hyperedges, so a
  printed or plotted snapshot says "cases" and "arbitrator" instead of
  "edge" and "node".
* `hg_edge_centrality()` and `hg_motifs()` on a temporal hypergraph report
  the snapshot `time` as a date for a calendar hypergraph and as the number
  on the clock otherwise; it used to be a character label of the day
  offset (`"3104"`), also when `at` was given as dates.
* `plot()` on a hypergraph with repeated member sets accepts `node_groups`:
  grouping the nodes now draws every hyperedge, as `color_by` does, instead
  of failing on the node sizes of the distinct-set view.
* `icsid_tribunals`: `respondent` and `subject` are trimmed of the spaces
  the source archive pads some values with. "Argentine Republic " held 12
  of Argentina's 47 cases, so a filter on the clean name missed them; there
  are now 145 distinct respondents instead of 201. No other column changes.
* `hg_subset(where =)` raises `hypergraphs_bad_input` for a value that no
  hyperedge takes, naming the closest values, instead of returning an empty
  hypergraph.

# hypergraphs 0.6.3

* `plot()` on a hypergraph gains `node_groups`: each node is coloured and
  shaped by its group, from a community fit (`hg_communities()`,
  `hg_mmsbm()`, `hg_topics()`) or any node/label table, so communities can
  be plotted over the hyperedges. The hulls turn grey unless `color_by` is
  given.
* `plot()` with `color_by = "size"` lists only the sizes of the hyperedges
  it draws, and the curves of a distribution over calendar times are named
  by date.
* A clock column of a hyperedge (`start`, `end`, `time`, `session`, ...) is
  no longer chosen as the count that titles and colours the hulls; the
  decision figures of a temporal snapshot are readable again.
* `hg_get()` on a hypergraph has a `weight` column only for a hypergraph of
  windows; other hypergraphs no longer print a column of `NA`.
* An unknown `what` in any accessor or verb raises `hypergraphs_bad_input`
  naming the available tables.
* `hg_agreement()` reads `hg_mmsbm()` and `hg_topics()` fits and text
  hypergraphs directly, and a table keyed by `doc` needs no `node =`.
* `hg_get()` on a community fit gains `converged`, `sort_by` and `top`.
* New article "Legal hypergraphs" on the tribunals of ICSID and the
  citations of the German Federal Constitutional Court, with hypergraph
  modularity.

* New dataset `forum_posts`: 360 simulated forum posts by 30 students whose
  sequence of topics has second-order memory, generated by
  `data-raw/forum_posts.R`. It is the example corpus for `hg_sequences()`.

* `read_hif()` and `write_hif()` are removed, and with them the `jsonlite`
  suggestion. A hypergraph is exchanged with other packages as a `netobject`
  or a `cograph_network`. `hg_get()` on a hypergraph no longer offers
  `what = "node_data"` or `"incidence_data"`, which only a read HIF file
  filled.

# hypergraphs 0.6.2

* `text_hypergraph()` gains `separator`: a delimited field such as the
  author keywords of a bibliographic export is read as whole-phrase terms,
  each a hyperedge over the papers that carry it, the keyword incidence of
  co-word analysis (Callon et al. 1983). Terms are lowercased and trimmed
  of the punctuation exports leave at their edges.
* New `hg_dictionary()` labels each node by the dictionary category whose
  terms it carries most (Grimmer & Stewart 2013); the result is the
  `labels` input of `hg_classify()` and `hg_hypergat()`.
* Documents without a label are `split = "unlabelled"` in an
  `hg_classification`, so `hg_get(fit, split = "unlabelled")` reads the
  classifier's proposals.
* `hg_get()` on a text hypergraph gains `node` (documents by id or by a
  table with a `node` column), `sort_by` and `top` (vocabulary).

* `hg_classify()` and `hg_hypergat()` gain `holdout`: a share of the labels,
  drawn within each class with `seed`, is hidden, predicted and scored. The
  result is an `hg_classification` that prints the held-out accuracy and
  balanced accuracy (Brodersen et al. 2010). `hg_get()` reads its
  `"predictions"`, `"accuracy"`, `"classes"` and `"confusion"` tables, and
  filters the predictions with `split`, `correct`, `node`, `sort_by` and
  `top` (per class). `plot()` draws the confusion table. A `hg_hypergat()`
  result retains its trained network. `hg_get()` computes attention on
  request, and `plot(type = "attention")` shows its diagnostic weights.
  `plot(type = "hyperedges")` plots a document's word hypergraph.
* `labels` of `hg_classify()` and `hg_hypergat()`, `clusters` of
  `hg_keywords()`, `topic_network()` and `hg_topic_sizes()`, and `group` of
  `hg_get(what = "prevalence")` accept the name of a document column
  (`labels = "subject"`).
* `plot()` on an `hg_topics()` model gains `type = "prevalence"` with
  `group`, the prevalence of every topic in each group.

* `hg_hypergat()` now returns a reusable fitted classifier in every run.
  `predict()` classifies new documents using its trained network, frozen
  vocabulary and topic keywords; optional observed labels evaluate these
  predictions without updating the model. `hg_get()` reads predictions,
  evaluation, document text, training history and optional attention
  diagnostics. `plot()` shows confusion or training loss. The vocabulary
  is estimated from known labels only. Unscorable held-out and new labelled
  documents count as errors, and prediction preserves all new input rows.
* HyperGAT attention diagnostics include uniform normalization baselines
  and the fraction of words exclusive to each edge. Attention plots show
  the baseline alongside the learned weights. Attention is computed only
  when requested and is described as an internal weight, not a prediction
  explanation. Legacy `what` extraction remains available with a warning.

* `hg_subset()` gains `component = "largest"`, which keeps the largest
  connected component, the input that label spreading and spectral
  clustering need. A subset text hypergraph now carries a matching text
  layer, so `hg_keywords()` and the other text verbs work on it.

* `stop_words_en()` gains `type = "snowball"`, the 174-word Snowball English
  stop list with pronouns and contractions, for chat and other informal text.
* `clean_text()` repairs emoji and other characters garbled by a
  UTF-8-as-Windows-1252 export.
* Labels or groups that name documents `text_hypergraph()` or
  `hg_hypergat()` dropped as empty are set aside with a warning in
  `hg_keywords()`, `hg_topic_sizes()`, `topic_network()`, `hg_classify()`,
  `hg_neural()` and `hg_hypergat()`, so the corpus table can be passed back
  whole; an id that never was a document is still an error.

* `hg_get()` on an `hg_topics()` model gains `what = "prevalence"` with
  `group =`: the mean share of each topic within each group of documents,
  the topic prevalence by covariate of Roberts et al. (2014).

* `clean_text()` gains `boilerplate = TRUE`, which removes publisher names,
  company suffixes and the phrases of licence and rights notices wherever
  they occur in a text, so a classifier of bibliographic abstracts does not
  learn the publisher and the year from the notice.

* `clean_text()` removes the word "Copyright" together with the notice that
  follows it ("Copyright © 2020 Elsevier Ltd." left "Copyright" behind).

* `plot()` on a `simplicial()` result gains `type =`. `"simplices"` (the
  default) draws the maximal simplices as before; `"summary"` draws the
  face counts, the Betti numbers and the simplicial degree, the summary
  figure of the complex.

* **Input formats with their own vocabulary.** Membership data name a
  `node` and a `hyperedge` (or `from` and `to`); event data name an
  `action` with its `session` and `actor`, as in the memory family;
  `group` is always the comparison variable, as in Nestimate. In
  `group_hypergraph()` and `temporal_hypergraph()`, `group =` (the
  hyperedge column) is now `hyperedge =` and `by =` is now `group =`.
  `hypergraph()` recognises the format from the arguments: `action` with
  `session` (or `actor`) gives one hyperedge per session, `window` gives
  windows, `node` with `hyperedge` gives membership hyperedges.
* **`pairwise_network()` builds the pairwise network of a hypergraph.** It
  returns a network object (`net_hg_pairwise`, also `netobject` and
  `cograph_network`) for every projection, chosen with `type = "clique"`
  (default), `"association"` or `"citation"`; `plot()` plots it through
  cograph, `hg_get()` returns its edges (`what = "nodes"` its nodes), and
  `hypergraph()` reads it as a network. It replaces `hg_project()` and
  `hg_clique_expansion()`, which are removed.
* **The topic network is `topic_network()`**, a constructor like
  `pairwise_network()`.
  It replaces `hg_network()`, its name in 0.6.1, and `hg_relations()`,
  which are removed: `hg_relations(hg, clusters)` is
  `hg_topic_network(hg, clusters = clusters)`.
* **Constructors are bare nouns.** `random_hypergraph(type = "uniform" |
  "regular" | "gnp" | "sbm")` replaces `hg_sample_uniform()`,
  `hg_sample_regular()`, `hg_sample_gnp()` and `hg_sample_sbm()`, with the
  same models and the same seeded draws; its arguments are given by name.
  `read_hif()` and `write_hif()` replace `hg_read_hif()` and
  `hg_write_hif()`.
* **`group_hypergraph(min_share =)`** keeps the frequent sets that reach a
  minimum support, the share of a group's sets that contain them (Agrawal &
  Srikant 1994); without `top`, every such set is kept. It applies to data
  frames, clustered sequences and topic models, and `print()` states the
  rule.
* `plot()` of a hypergraph packs disconnected pieces into one frame by
  default (`pieces = "packed"`), so each piece takes room in proportion to its
  size; `pieces = "row"` gives the former side-by-side frames.
* `plot()` of a hypergraph whose hyperedges repeat (the trials of an event
  log) draws its distinct sets, each coloured by its number of copies; a
  hypergraph of counted sets with one group is drawn as that group; a window
  hypergraph is coloured by its window counts. `plot()` of a
  `pairwise_network()` uses a circle layout. None of these needs an argument.
* `plot()` of a hypergraph takes a numeric hyperedge attribute as its count
  only when the attribute varies between hyperedges; a constant one, such as
  the session number of the runs of one session, no longer colours and titles
  the hyperedges. `color_by = "weight"` colours the hyperedges of a
  `window_hypergraph()` by their window counts.
* **`tutoring_events` replaces `debug_events`.** It holds every step of the
  tutoring data, 84,356 events from 13,309 problem steps (20,626 trials),
  with each event renamed to a word of the same meaning and similar events
  merged (19 events become 15). `debug_events` held a random 4,000 steps
  under a coding-assistant vocabulary. Its `outcome` column marks a step
  `completed` (with a correct answer) or `stopped` (without one).
* `hg_measures()` builds the hyperedge-by-hyperedge overlap matrices only for
  `what = "overlap"`. The node table and the summary of a hypergraph with
  20,000 hyperedges took minutes or exhausted memory; they take under a
  second.
* `hg_subset()` keeps the window counts of a `window_hypergraph()` aligned
  with the hyperedges it keeps; printing or plotting such a subset failed.
* `hg_agreement()` accepts a fit of `hg_communities()` and compares its
  medoid partition.
* New vignette, *Hypergraphs*: observed groups, frequent sets, windows,
  projection, node measures, centrality, communities and null models on
  `debug_events`.

# hypergraphs 0.6.1

* **The package is renamed from hypernets to hypergraphs.** Every verb keeps
  its name; condition classes and result classes change their prefix from
  `hypernets_` to `hypergraphs_` (`hypergraphs_bad_input`,
  `hypergraphs_result`, ...), so code that catches a condition by class
  changes with them.
* **`hypergraph()` is the main constructor.** It reads any input and calls
  the constructor of that kind of data, returning its result unchanged: a
  data frame of groups goes to `group_hypergraph()` (with `by` or `top`, the
  frequent sets), with a clock to `temporal_hypergraph()`, with `window`,
  `step` or `action` to `window_hypergraph()`; a list of sequences to
  `window_hypergraph()`; a network (matrix, sparse matrix, `netobject`,
  `cograph_network`) to `network_hypergraph()`; a topic model or a clustering
  of sequences to `group_hypergraph()`. The specific constructors remain.

* **`group_hypergraph(by =)` counts frequent sets within groups.** With
  `by`, each value of `group` (a trial, a session, a basket) is one set of
  `actor` values, and the `top` most frequent sets within each value of `by`
  (an outcome, a cluster) become hyperedges, read with
  `hg_get(x, what = "sets")` and drawn per value with `plot(x, group = )`,
  as for clustered sequences. The topic-combination route names its grouping
  column `by` as well. Giving `top` without `by` counts the sets over all
  groups. New dataset `debug_events`: 4,000 sessions with a coding assistant,
  one row per event, with runs (a stretch of a session between two verdicts)
  and a quick or slow group; real event sequences of another domain with
  every event renamed (`data-raw/debug_events.R`).
* `plot()` of a hypergraph colours `color_by = "size"` discretely, one
  Okabe-Ito colour per size; it was a continuous scale.
* `network_hypergraph(type = "vr")` raises `hypergraphs_bad_input` and points to
  `simplicial(type = "vr")`; it raised a plain error that said the package
  builds no Vietoris-Rips complex.

* **Communities of memory networks can be drawn over the transition network.**
  `plot()` of a `hg_communities()` result can draw the network as a transition
  network analysis plot (TNA styling through
  `cograph::overlay_communities()`): arrows carry transition probabilities,
  node area follows flow, and every community is a blob around its nodes.
  New `type = "network"` draws the network of states, in which a shared state
  lies in two blobs; `type = "states"` now draws the memory nodes the same
  way, laid out by the walk's transitions between them, in place of
  cograph's group-ring layout. The default pebble view of the community
  hypergraph is unchanged, and the labels say "state" and "memory node" as the
  tables do.

* `hg_network()` replaces `hg_relations()` (kept as a deprecated alias
  that warns with `hypergraphs_deprecated` and returns the identical
  result). Besides a partition (`clusters =`), it takes a topic model
  (`topics =`). Without `threshold` two topics are linked by the
  correlation of their shares (stm's simple `topicCorr()`, equal on the
  local oracle); with `threshold` the weight is the number of documents in
  which both topics reach that share (Abuhay et al. 2017; Cassi et al.
  2017), normalisable by the same similarity measures.
* `group_hypergraph(topic_model, threshold =)` builds the hypergraph of
  topic combinations: each document is the set of its topics at the
  threshold, and the most frequent sets become hyperedges, counted as for
  clustered sequences (`hg_get(x, what = "sets")`). `by =` counts them
  within a column of the documents table (a year, an author); `min_size`
  and `top = Inf` keep the multi-topic sets and all of them. Hypergraphs of
  clustered sequences now print their source.
* `hg_sequences(topics =)` builds sequences from each document's main
  topic, grouped by any column of the documents table.

* `hg_markov_stability()` reads only the possible transitions of a matrix
  when it checks irreducibility, so an unnormalised matrix warns once with
  `normalize = TRUE` (it warned twice) and is refused without a warning with
  `normalize = FALSE` (it warned that it was normalising).

# hypernets 0.6.0

* **One naming grammar and one reader.** Constructors (and `markov_order()`, `memory()`) are bare nouns,
  every other verb is `hg_*()`, and every table of every object is read with
  `hg_get(x, what = )`. No `build_*()`, no `hon_*()` and no `as.data.frame()`
  method remains; no export shares a name with an imported package.
  Every number is unchanged: the memory and simplicial verbs return the
  estimator's objects as before, and each old accessor table is
  `identical()` to its `hg_get()` replacement (58 tables checked against
  the previous build; the one difference is `hg_degree()`'s reset row
  names).

  | Before | After |
  |---|---|
  | `as.data.frame(x, what = )` (every class) | `hg_get(x, what = )` |
  | `build_hon()` | `hon()` |
  | `build_honem()` | `honem()` |
  | `build_mogen()` | `mogen()` |
  | `build_hypa()` | `hypa()` on sequences (on a hypergraph: the pair HYPA, as before) |
  | `build_simplicial()` | `simplicial()` (now verifies a clique complex; see below) |
  | `markov_order_test()` | `markov_order()` |
  | `path_dependence()` | `memory()` |
  | `markov_stability()` | `hg_markov_stability()` |
  | `persistent_homology()` | `hg_homology()` |
  | `persistence_landscape()` | `hg_landscape()` |
  | `bottleneck_distance()` | `hg_bottleneck()` |
  | `wasserstein_distance()` | `hg_wasserstein()` |
  | `betti_numbers()` | `hg_betti()` |
  | `euler_characteristic()` | `hg_euler()` |
  | `q_analysis()` | `hg_qanalysis()` |
  | `simplicial_degree()` | `hg_degree()`, or `hg_get(sc, what = "degree")` |
  | `verify_simplicial()` | folded into `simplicial(verify = TRUE)` |
  | `bootstrap_hon()` | `hg_bootstrap()` |
  | `compare_hon()` | `hg_compare()` |
  | `hon_centrality()` | `hg_centrality()` on a memory network |
  | `hon_communities()` | `hg_communities()` on a memory network |
  | `mogen_transitions(mg, order =)` | `hg_get(mg, what = "transitions", order =)` |
  | `path_counts(data, k =)` | `hg_get(mogen(data), what = "paths", k =)` |
  | `pathways(x)` | `hg_get(x, what = "pathways")` (a one-column table, `pathway`) |

  - **`hg_get(x, what = NULL, ...)`** is an S3 generic with a method for
    every result class, including the memory and simplicial classes;
    `what = NULL` is the primary table, and filters, `sort_by` and `top`
    keep their names. An object without a method raises
    `hypernets_bad_input` naming its class. Print methods now point at
    `hg_get()`.
  - **One verb per idea.** `hg_centrality()`, `hg_communities()` and
    `hypa()` dispatch on the class of their input: a memory network
    (`net_hon`, or sequences for `hypa()`) or a hypergraph (`net_hg`).
    Each method keeps its arguments and defaults; an argument that only the
    other method takes raises `hypernets_bad_input` instead of being
    ignored. The first argument of the three is now `x` (was `hg` / `hon`).
  - **Every sequence-taking verb reads every input form**: a long event
    table (`action =`, `actor =`, `time =`, and new `session =`), a wide
    data frame, a list of sequences, or a model object (netobject,
    netobject_group, tna, cograph_network). `hon(long, action =, actor =,
    time =)` is back and is `identical()` to building the relative
    transition network with an infinite session gap and passing it on. A
    position column is passed as `time`. A long table passed without
    `action =` raises `hypernets_long_format` again, now recognising the
    usual column names by role (`code`/`state`/`event`, `session_id`/`user`,
    `timestamp`/`order`, ...).
  - **`hg_get(mg, what = "paths", k =)`** reads k-state path counts off the
    fitted order-(k - 1) layer of a `mogen()` model (`k` up to the highest
    fitted order plus one); it equals the old raw-data counter on sequences
    without internal gaps. `hg_get(x, what = "pathways")` works on `hon()`,
    `mogen()` and sequence `hypa()` fits.
  - **`simplicial()` checks a clique complex on construction**: its
    simplices against igraph's cliques of the same thresholded graph (when
    igraph is installed) and its Euler characteristic against the
    alternating Betti sum. A failure raises the warning
    `hypernets_simplicial_unverified`; `verify = FALSE` skips the check.
  - Dropped: `hypa(k =)`, the deprecated alias of `order`, now raises
    `hypernets_bad_input`.

* **One vocabulary for table columns, and prints that show the table.**
  Every table `hg_get()` returns names the same quantity the same way:
  `count`, `expected`, `z`, `p_value`, `p_adj`, `significant`,
  `log_likelihood`, `df`, `aic`, `bic`, `community`, `run`, `dimension`,
  `node`, `members`. Renamed columns (old -> new):
  - `markov_order()`: `loglik` -> `log_likelihood`, `AIC`/`BIC` ->
    `aic`/`bic`, `p_permutation` -> `p_value`. `mogen()`: `dof` -> `df`,
    `layer_dof` -> `layer_df`. `hg_markov_stability()`: `stationary_prob`
    -> `stationary`. `memory()`: `n` -> `count`, `H_order1`, `H_orderk`,
    `H_drop` -> `entropy_first_order`, `entropy_order_k`, `entropy_drop`,
    `KL` -> `kl` (and `sort_by = "kl"`), `top_o1`/`top_ok` ->
    `top_first_order`/`top_order_k`. `honem()` variance: `dim` ->
    `dimension`. `hon()` nodes: `id`, `node` (the duplicated `label` and
    `name` are gone).
  - `hypa()`: `observed` -> `count`, `p_tail` -> `p_adj`, `anomaly` ->
    `direction`, `p_adjusted_under`/`p_adjusted_over` ->
    `p_adj_under`/`p_adj_over`. `hg_compare()`: `prob_<group>` ->
    `probability_<group>`.
  - `hg_communities()` on a memory network: the memory node is `node`, its
    physical state `state`, the cluster `community` (was `state`,
    `physical`, `module`); `n_modules` -> `n_communities`, `trial` -> `run`,
    `n_states`/`n_physical` -> `n_nodes`/`n_states`; the filter argument is
    `community =` and a given `partition` data frame has `node` and
    `community`. The dataset `ring_communities` has `node` for `state`.
  - Hypergraphs: the edge table's `states` -> `members`; temporal
    memberships `member` -> `node` (`edge` stays, as in HIF);
    `hg_mmsbm()` `u` -> `membership_weight`, `w` -> `affinity`, restarts
    `start`/`loglik`/`ari_best` -> `run`/`log_likelihood`/`ari_to_best`;
    `hg_motifs()` `observed`/`null_mean`/`n` -> `count`/`expected`/`n_null`;
    `hg_compare_communities()` `largest`/`second` ->
    `largest_size`/`second_size`, and `what = "matrix"` is gone (the
    `similarity` table holds the same values).
  - Text hypergraphs: `n` -> `count` in the weights and vocabulary tables.
  - Simplicial complexes: `id`/`dim` -> `simplex`/`dimension`, the filter
    argument `dim =` -> `dimension =`; the validation table `set` ->
    `members`, `validated` -> `significant`.

  Every result prints a header line with its main settings, then the first
  rows of its default table (`print(x, n = )` for more). The results of
  `hon()`, `honem()`, `mogen()`, `markov_order()`, `memory()`,
  `hg_markov_stability()`, `hg_homology()`, `hg_landscape()` and
  `hg_qanalysis()` carry the class `hypernets_result` in front of the
  estimator's classes for this; plot() and the estimator's own verbs are
  the estimator's.

  `summary()` on a result returns every table of the result as a list of
  data frames named after the `what` values of `hg_get()`, so
  `summary(x)$validation` is the data frame `hg_get(x, what =
  "validation")` returns; overall figures are added as further tables
  (`overall` for `memory()`, `markov_order()` and `hg_compare()`,
  `by_order` for `hg_bootstrap()` and `hg_compare()`, `communities` for
  `hg_mmsbm()`). The result itself is not changed. This replaces the
  printed reports of the estimators' summaries and the single tables the
  bootstrap, comparison, community and mixed-membership summaries
  returned.

* **Sequences are built exactly as in the tna family.** Every verb that
  reads a long event table (`hon()`, `mogen()`, `hypa()`, `markov_order()`,
  `memory()`, `hg_markov_stability()`, `hg_bootstrap()`, `hg_compare()`
  through the group model, `simplicial(type = "window")`,
  `window_hypergraph()`) builds the sequences that
  `Nestimate::build_network(method = "relative")` builds from the same
  arguments (`identical()`, tested on six call patterns of `human_long`):
  - new `time_threshold = 900` and `timezone = "UTC"`: with `time`, a gap
    of more than 900 seconds starts a new sequence; `time_threshold =
    FALSE` keeps one sequence per actor or session. On `human_long` ordered
    by `timestamp` this gives 526 sequences instead of 429.
  - a column named `action`, `time`, `session` or `session_id` is used when
    its argument is `NULL`; `session = FALSE` switches session detection
    off (with `actor = "project"`, 34 sequences instead of the 429 sessions
    detected).
  - `actor` and `session` may name several columns.
  - a missing actor or session raises `hypernets_bad_input` (rows were
    dropped before); without `actor`, the message
    `hypernets_single_sequence` says all events form one sequence.
  - `window_hypergraph()` gains `session`, `time_threshold` and `timezone`.

* **`hg_markov_stability()` refuses a chain that is not irreducible.** It
  used to return return and passage times of the order of 10^15 steps with a
  plain warning. A chain with a transient or absorbing state, several closed
  classes, or a state with no outgoing transition now raises
  `hypernets_not_ergodic` naming those states. `plot()` of a stability
  result draws one panel of bars per measure, every state included (the
  estimator's plot dropped states beyond its eight colours), and
  `plot(x, what = "passage_time")` the first-passage heatmap through cograph.

* **`hg_topics()`: a mixed-membership topic model.** The document-word
  counts of a hypergraph are factorized by the multiplicative updates that
  minimise the Kullback-Leibler divergence (Lee & Seung 1999, 2001), which
  is probabilistic latent semantic analysis (Hofmann 1999; Gaussier &
  Goutte 2005): every document is a mixture of topics and every topic a
  distribution over words. The fit keeps the best of `nstart` seeded starts
  (`parallel = TRUE` gives the same result) and reports each topic's
  agreement across the starts by the matched average Jaccard of its top
  words (Greene, O'Callaghan & Cunningham 2014), topics matched by the
  Hungarian method. `hg_get(fit, what = "topics" | "shares" | "words" |
  "documents" | "restarts")`, with `print`, `summary` and `plot`. The
  updates agree with scikit-learn's KL solver to 1e-11 after 200 steps
  (local oracle).
  - `hg_topic_quality(hg, topics = fit)` scores a topic model: UMass or NPMI
    coherence of its most probable words and FREX from its word
    distributions, both identical to stm on the same beta (local oracle).
  - `hg_cluster(algorithm = "symnmf", what = "membership")` returns the
    graded memberships of the symmetric factorization (Kuang, Ding & Park
    2012), each node's row normalised to sum to one.
  - New dataset `covid_sample`: a seeded simple random sample of 1,000
    abstracts from the COVID-19 education corpus of the sbert package
    (`data-raw/covid_sample.R`), and a new vignette, "Mixed-membership
    topics of the COVID-19 education literature", built on it. The
    existing vignette and `covid_abstracts` are unchanged.

* **Group models of memory networks.** `hon(data, ..., group =)` takes a
  column name or one label per sequence and returns one memory network per
  group (`net_hon_group`); `hg_get()` stacks their tables with a `group`
  column. `hg_compare()` takes the group model (`groups =` picks two when
  there are more) instead of two data sets and repeated settings, and
  `hg_bootstrap()` on a group model bootstraps every group. Sequence
  `hypa()` gains `type = c("all", "over", "under")`, `order_by = c("sig",
  "ratio", "freq", "path")` and `n` for its printed view, and
  `hg_get(fit, what = "over" / "under", order_by =, top =)` sorts each
  direction by its own tail.

* **Simplicial complexes from sequences.**
  - `simplicial(type = "window", window =, min_count =)` builds the
    co-occurrence complex of windows of consecutive actions (Salnikov et al.
    2018), from the same sequence input as the memory verbs. A set of
    actions seen in at least `min_count` windows is a simplex with all its
    faces.
  - `simplicial(type = "window", validate = TRUE)` needs no count
    threshold: a set of actions becomes a simplex when it occurs in more
    windows than a null model predicts, with Benjamini-Hochberg control over
    all possible sets of each size (`alpha = 0.05`). The default
    `null = "swap"` is swap randomization (Gionis et al. 2007): every window
    keeps its number of actions and every action its number of windows.
    `n_null = NULL` takes enough shuffles for a single set of any size to
    pass (at least 999); fewer raise `hypernets_low_resolution`. On
    `human_long` (window 3, 1680 shuffles) 30 sets pass under each of five
    seeds, 2 more under some, and the Betti numbers are (1, 3, 0) under all
    five. `null = "hypergeometric"` is the statistically validated
    hypergraph of Musciotto, Battiston & Mantegna (2021), with the exact
    p-value of their Eq. 3; its null lets an action fall into windows
    independently of the others, so with few actions it expects more
    co-occurrence than the windows hold (1.7 times for pairs and 2.9 times
    for triples on `human_long`, where no set passes, also when the
    sessions are pooled by project or into one sequence). It suits sparse
    co-occurrence, each set expected in less than one window, and warns
    `hypernets_dense_cooccurrence` otherwise.
    `hg_get(sc, what = "validation")` returns the tests.
  - `hg_homology()` of a window complex computes persistence over the count
    filtration (Petri et al. 2013): a simplex enters at the largest number
    of windows of any set containing it, and the Betti curve at count `t`
    is the Betti numbers of the complex with `min_count = t`.
  - `plot()` of a complex draws its maximal simplices, all of them, the
    most significant (validated), most frequent (window) or closest
    (Vietoris-Rips) first, in Okabe-Ito colours; `dismantled = TRUE` draws
    one panel per simplex and `top =` limits the number. It used to hand the
    complex to cograph, which drew faces as well and only the first ten.
  - `simplicial(type = "clique")` takes sequences with `actor`, `action`
    and `time` and builds their relative transition network.
  - `simplicial(type = "vr")` and `hg_homology(type = "vr")` take points (a
    non-square table or matrix, a data frame of coordinates, or a `dist`
    object) as well as distances.
  - `hg_get(sc, what = "betti")` tabulates the Betti numbers; `plot()` of a
    persistence landscape draws the landscapes that are not zero, with one
    colour and line type each.
  - `simplicial(type = "clique", direction = "both")` joins a pair only
    when both directed weights reach `threshold`, and then reproduces the
    cliques of `tna::cliques()` exactly (24 of 24 size and threshold
    combinations on two data sets). `"either"` (default) is the previous
    rule. The clique complex now cites Giusti, Ghrist & Bassett (2016).

* **The memory and simplicial families are Nestimate's; hypernets imports
  it.** hypernets now `Imports: Nestimate (>= 0.8.5)` and re-exports
  Nestimate's `build_hon()`, `build_honem()`, `build_hypa()`, `build_mogen()`,
  `markov_order_test()`, `markov_stability()`, `path_dependence()`,
  `mogen_transitions()`, `path_counts()`, `pathways()`, `build_simplicial()`,
  `persistent_homology()`, `q_analysis()`, `betti_numbers()`,
  `euler_characteristic()`, `simplicial_degree()`, `verify_simplicial()`,
  `bottleneck_distance()` and `persistence_landscape()` instead of carrying
  copies (the copies matched CRAN Nestimate 0.8.5 on 22 of 23 calls at
  tolerance 0; `markov_stability()` differed only by Nestimate's rounding).
  Nestimate does not depend on hypernets. hypernets keeps `as.data.frame()`
  for Nestimate's classes (Nestimate has none) and its own verbs on top of
  them: `bootstrap_hon()`, `compare_hon()`, `hon_centrality()`,
  `hon_communities()`, `wasserstein_distance()`. Loading both packages now
  masks nothing and overwrites no S3 method. Result classes are Nestimate's
  (`simplicial_complex`, `persistent_homology`, `q_analysis`,
  `persistence_landscape`; the 0.5.x `net_*` simplicial names are gone).
  Lost with the hypernets copies, because Nestimate lacks them:
  - `build_hon(action =, actor =, time =)` long input -- build the network
    with `Nestimate::build_network(data, method = "relative", actor =,
    action =, time =, time_threshold = Inf)` and pass it to `build_hon()`
    (identical matrix); `bootstrap_hon()` / `compare_hon()` still take the
    long form directly.
  - the long-format guard on `build_hon()`, `build_mogen()`, `build_hypa()`,
    `markov_order_test()` (a long table passed bare is read as wide).
  - `top =` on `mogen_transitions()` and `simplicial_degree()`
    (`as.data.frame(x, top =)` still truncates); validation of
    `path_counts(top =)`.
  - hypernets' `markov_stability()` rewrite: the reducible-chain error
    (`hypernets_not_ergodic`), classed input errors, list-of-sequences
    input, unrounded `$stability`, and its landscape / states / passage /
    network plots.

* **`hg_mmsbm()`: probabilistic (mixed) membership.** Hy-MMSBM (Ruggeri,
  Contisciani, Battiston & De Bacco 2023): every node (a document, in a text
  hypergraph) gets a probability of belonging to each of `k` communities,
  fitted by the authors' EM with `nstart` restarts; hyperedges of any size
  are fitted exactly (no truncation unless `max_size` is set).
  `as.data.frame(fit, what = "membership" / "nodes" / "affinity" /
  "restarts")`; the restarts table reports each start's log-likelihood,
  convergence and ARI with the kept partition, so instability is visible.
  Local oracle against the authors' code (commit 6a12077): one EM step to
  4.8e-15, log-likelihood to 1.1e-15, full fits to 2.1e-12. Two deliberate
  deviations: restarts are compared by the paper's Eq. 5 log-likelihood (the
  authors' CLI uses an inconsistent C = 1 form), and the default stopping
  rule is on normalised memberships (`criterion = "membership"`), because
  the authors' parameter-change rule never fires under their default priors
  (the likelihood is invariant to u -> cu, w -> w/c^2); theirs is
  `criterion = "parameters"`. New warning class `hypernets_isolated_nodes`.

* **Hypergraph clustering, assortativity and Katz centrality.**
  `hg_transitivity()` (local clustering: projection, Watts & Strogatz 1998;
  extra overlap, Zhou & Nakhleh 2011 / Klimm et al. 2021; two-node
  union / min / max, Latapy et al. 2008), `hg_assortativity()` (Chodrow
  2020: uniform, top-2, top-bottom; rank or degree scale),
  `hg_degree_correlation()` (Lotito et al. 2023), and
  `hg_centrality(type = "katz", alpha =)` (Katz 1953 on the
  Estrada & Rodriguez-Velazquez 2006 hypergraph adjacency). Local oracles:
  XGI 0.10.2 (to 2.2e-16), hypergraphx 1.8.0, igraph `alpha_centrality`.
  Undefined coefficients are `NA` (XGI reports 0); XGI's exact uniform
  assortativity weights hyperedges by m(m-1) and so differs from Chodrow's
  definition on mixed sizes (asserted). Estrada & Rodriguez-Velazquez's
  global C2(H) is not included (no oracle; see `workinprogress/`).

* **`hg_modularity()` and `hg_communities(type = "irmm")`.** Hypergraph
  modularity of a partition (Kaminski et al. 2019; `type = "linear"`,
  `"majority"`, `"strict"` as in HyperNetX) as a tidy score or per-community
  table, and IRMM community detection (Kumar et al. 2020) inside the existing
  `hg_communities()` ensemble (seeded runs, AMI medoid, comparison with
  Infomap through `hg_compare_communities()`; `as.data.frame(fit, what =
  "weights")` gives the reweighted hyperedges). Local oracle: HyperNetX
  2.4.3, 90 scores to 5.6e-16; IRMM passes replayed step by step to
  4.4e-16. (HyperNetX truncates fractional weights in its degree tax; our
  value equals HNX on the same weights rescaled to integers.)

* **`hg_write_hif()` / `hg_read_hif()`** read and write the Hypergraph
  Interchange Format (Coll et al. 2025): schema-valid (official JSON
  schema), doubles lossless, `window_counts` as edge weights, node / edge /
  incidence attributes kept (`as.data.frame(hg, what = "node_data" /
  "incidence_data")`); round-tripped through XGI 0.10.2 and HyperNetX
  2.4.3. Directed HIF is refused.

* **Fixes and small API changes found while writing the 0.6.0 documents.**
  - `hg_pagerank()` ignored a window hypergraph's `window_counts` and used
    the SD+1 heuristic, so it disagreed with `hg_centrality(type =
    "pagerank")` and the Laplacian walk (by 0.018 on `ring_sequences`,
    window 3). It now uses the package-wide default: explicit
    `edge_weights`, else `window_counts`, else the heuristic. The sparse
    walk operators follow the same rule (they matter for
    `hg_read_hif(sparse = TRUE)` files with edge weights).
  - `hg_write_hif()` wrote missing attribute values as `null`; they are now
    omitted, so read -> write is a fixed point for files where only some
    records carry an attribute. `as.data.frame(hg, what = "edge_data")` is
    new, and `print()` names HIF as the source of a read hypergraph.
  - `hg_topic_quality(words =)` also takes `hg_keywords()` output (its
    `cluster` column is the topic) and the character matrix
    `topicmodels::terms()` returns, as they are.
  - `as.data.frame(<hg_communities>, what = "ami" / "ari" / "nmi")` returns
    one row per distinct pair of runs (`run_a`, `run_b`, value) instead of
    the full matrix with `Var1` / `Var2` and the diagonal.
  - `hg_mmsbm()` normalised memberships that EM had driven to zero or to
    subnormal remnants, reporting exact 1/3, 1/2 ... "mixtures" for them,
    and left exact-zero rows `NA` without a warning. A node whose row total
    is below working precision relative to the largest row now has no
    membership (`NA`), is counted in `print()`, and raises
    `hypernets_collapsed_membership`. Subnormal affinities are returned as 0.
  - `hg_keywords()` and `hg_agreement()` accept a `community` column (the
    `hg_communities()` medoid, `hg_mmsbm()` node table). `hg_agreement()`
    left NA-labelled nodes in `n` and `agreement` while the indices dropped
    them; they are now dropped for every column, with a
    `hypernets_missing_labels` warning.
  - `?hg_assortativity` states that `"top_2"` and `"top_bottom"` are
    positive on random hypergraphs (order statistics), so 0 is not their
    null value.

* **Bug fix: the neural verbs lost half of every symmetric matrix.**
  `.thg_torch_sparse()` coerced to `"TsparseMatrix"` and read its triplets;
  Matrix stores a symmetric matrix as `dsTMatrix` (upper triangle only) and a
  unit-triangular one without its diagonal, so those cells were dropped
  silently. `hg_hypergcn()` therefore propagated through an
  upper-triangular adjacency in all three methods (logits off by up to 0.83
  against the official HyperGCN code), and any symmetric `features` matrix
  given to `hg_hnhn()`, `hg_allset()` or `hg_neural()` was truncated. The
  conversion now expands to a general matrix first. Found by the new DHG /
  official-code oracle (`local_testing_and_equivalence/test-oracle-hypergcn-hnhn.R`),
  which now agrees to about 6e-8.

* **Topic measures match stm; stability resamples documents.**
  - `hg_topic_quality()` now defaults to UMass coherence (Mimno et al. 2011)
    and FREX exclusivity (Bischof & Airoldi 2012) computed exactly as
    `stm::semanticCoherence()` / `stm::exclusivity()` (local oracle: max
    difference 0 on 80 abstracts); `coherence = "npmi"` gives corpus NPMI
    (Bouma 2009; Lau et al. 2014; text2vec oracle to 1e-8). The previous
    within-cluster NPMI and mean share are `coherence = "npmi_cluster",
    exclusivity = "share", sort_by = "share"` (identical output). `words =`
    scores topic word lists from another model (e.g. topicmodels, stm) with
    the same measures. New columns `coherence_type`, `exclusivity_type`.
  - `hg_stability()` defaults to `resample = "subset"`: Hennig's (2007)
    subsampling stability, the per-cluster best-match Jaccard over `n_boot`
    node subsamples (fpc `clusterboot(bootmethod = "subset")` oracle, max
    difference 1.1e-16), plus the eigengap (von Luxburg 2007) per `k`. The
    old two-seed comparison, which never resampled the data and so measured
    only solver determinism, is `resample = "seeds"`.
  - `hg_cocluster()` (new): spectral co-clustering of a bipartite incidence,
    nodes and hyperedges together (Dhillon 2001; sklearn
    `SpectralCoclustering` oracle).
  - `hg_keywords()` on a word-node bag hypergraph (`nodes = "word"`) raised
    no error and returned document ids in its `word` column; it now raises
    `hypernets_bad_input`.

* **Hypergraph verbs are `hg_*()`; the class is `net_hg`.** hypernets will
  import Nestimate (ROADMAP Phase 0b), which keeps its own
  `build_hypergraph()`, `hypergraph_*()` verbs and `net_hypergraph` class, so
  both packages load in every session. Shared names would mask each other's
  verbs and overwrite each other's S3 methods (measured with both loaded: 26
  masked, 32 overwritten). hypernets was never released, so the renames
  below break no published code.

  | Was | Is |
  |---|---|
  | `build_hypergraph()` | `network_hypergraph()` (the `<source>_hypergraph()` rule) |
  | `clique_expansion()` | `hg_clique_expansion()` (Nestimate's `clique_expansion()` accepts only its own class) |
  | 27 `hypergraph_*()` aliases of `hg_*()` (`hypergraph_pagerank`, `hypergraph_motifs`, ...) | removed; use the `hg_*()` name |
  | `hypergraph_allset()`, `hypergraph_hypergcn()`, `hypergraph_hnhn()`, `hypergraph_snapshot(s)()`, `hypergraph_laplacian()`, `hypergraph_joint_cluster()` | `hg_allset()`, `hg_hypergcn()`, `hg_hnhn()`, `hg_snapshot(s)()`, `hg_laplacian()`, `hg_joint_cluster()` |
  | `hypergraph_alldeepsets()`, `hypergraph_allset_transformer()` | `hg_allset(model = "deepsets" / "transformer")` |
  | engines `hypergraph_centrality()`, `hypergraph_cluster()`, `hypergraph_measures()`, `hypergraph_transduction()` | internal; `hg_centrality()`, `hg_cluster()`, `hg_measures()`, `hg_classify()` are the verbs and now carry the engines' documentation and references |
  | class `net_hypergraph` (+ `_cluster`, `_transduction`, `_measures`, `_snapshots`) | `net_hg` (+ the same suffixes) |

  `hg_centrality()` gains the engine's `"pagerank"` and `"subhypergraph"`
  types with `damping` and `edge_weights`; `hg_classify()` gains
  `edge_weights`. Both are `identical()` to the engine (new tests). No
  computed value changes: the R8 transduction guard still gives 0.8451, and
  the Nestimate identity tests pass with a names-only normalizer.
  `tests/testthat/test-api-names.R` fails if an export or S3 registration
  reintroduces a Nestimate hypergraph name.

# hypernets 0.5.1

* **`group_hypergraph()` reads Nestimate clusterings directly.** A mixture
  Markov fit (`net_mmm`), a distance clustering (`net_clustering`) or the
  per-cluster networks built from either (`netobject_group`, with names from
  `rename_models()`) become a hypergraph of each cluster's most frequent
  state sets: every sequence reduces to the set of its distinct states and
  the `top` (default 8) most frequent sets of each cluster are the
  hyperedges, carrying `group`, `set` and `count` (frequent-itemset support
  counting; Agrawal & Srikant 1994). `states =` keeps a subset of states
  before counting. `as.data.frame(hg, what = "sets")` and
  `what = "state_counts"` return the tables; `plot(hg, group = "Cluster 1")`
  draws one cluster's sets, nodes sized by the cluster's sequences
  containing the state, sets named by what they add to the states all of
  them share. The objects are read by structure; Nestimate is not a
  dependency. On the Eventdata26 three-cluster mixture (750 / 1,784 / 1,301
  steps) the sets and counts are identical to the pipeline's
  `step_variations()` for every cluster (local parity test). The data.frame
  path is unchanged (`identical()` to a fixture frozen before the change).

* **`markov_order_test()` cites its sources.** New `@references` and a
  `@details` section tying each output column to its origin: the
  likelihood-ratio (`g2`, `df`, `p_asymptotic`) test of Anderson & Goodman
  (1957); the within-context permutation (`p_permutation`) as the
  margin-fixed conditional test of Agresti (1992), with the caveat that
  overlapping tuples from one trajectory are serially dependent and the
  exact test for a chain conditions on transition counts (Besag & Mondal
  2013); the multi-order log-likelihood of Scholtes (2017) and AIC/BIC order
  selection (Tong 1975; Katz 1981). Two implementation choices are stated as
  such: `df` counts observed categories per context, and layer parameters
  are counted on observed transitions. No computed value changes.

* **`betti_numbers()` and `persistent_homology()` document their
  coefficient fields.** `betti_numbers()` ranks oriented boundary matrices
  over the rationals (`qr()`); `persistent_homology()` reduces over Z/2.
  They differ when integral homology has torsion: on the 6-vertex real
  projective plane `betti_numbers()` gives (1, 0, 0) and
  `persistent_homology()` has essential classes (1, 1, 1). A test pins both
  (Hatcher 2002). No computed value changes.

* **HyperGAT benchmark compares like with like.** The `hg_hypergat()`
  benchmark rows (R8 0.9665, R52 0.9433) were run with `semantic = "none"`
  and are now set against the paper's "w/o semantic" ablation (Ding et al.
  2020, Table 4: R8 0.9714, R52 0.9415) rather than full HyperGAT
  (Table 2: 0.9797, 0.9498). The result CSV records `semantic`; the
  published table gains the ablation rows; `RESULTS.md`,
  `render_results.R` and the benchmarks article compute the comparison
  from them. No training was rerun.

* **`plot.net_hypergraph()` draws node overlays** as ordinary ggplot layers.
  `node_sizes =` (a `node`/`value` table or a vector named by node) draws each
  node as a circle in data units whose area follows the value (largest radius
  4% of the layout, none below 30% of it, so rare nodes stay visible);
  `direction =` (a `from`/`to`/`weight` table) adds a triangle pointing at the
  node that most often follows it, ties broken by name, none for a sink;
  `arrow_style = "inside" | "outside"`, `node_fill`, `arrow_fill`;
  `transitions =` draws the moves as curved arrows with width by `weight` and
  self-loops; labels move above each circle with a white halo, and the
  caption states what circle area (`size_title`) and triangles show.

* **`plot.net_hypergraph()` draws the event-data blob figures exactly.**
  Every figure the Eventdata26 pipeline drew with its `plot_blobs()` helper
  (19 calls across six documents) is now one `plot()` call on
  `group_hypergraph(members, actor = "state", group = "group")`, with the
  same built layer data (maximum absolute difference 0 in every layer;
  local-only parity test), and the call passes only the data:
  `plot(groups, node_sizes = node_sizes, notes = notes)`. A numeric
  hyperedge attribute that `group_hypergraph()` kept (such as `trials`,
  constant within each group) colours the pebbles, writes a title box with
  its count beside each one and names the unit, with no argument; of several
  numeric attributes the one that varies between hyperedges is read.
  `color_by`, `titles`, `unit` and `legend_title` override it. That look is
  the default for every hypergraph plot: haloed bold labels, legends below
  (2.5 cm keys for the colour bar, stacked when transitions add a third;
  discrete legends keep default keys, at most four to a row), margins
  40/130/30/130 pt, `alpha = 0.5`, `label_size = 4.2`, and `pieces = "row"`
  (each disconnected piece laid out on its own and set side by side;
  `pieces = "packed"` packs them). The earlier look is partly available
  through `alpha = 0.45`, `label_size = 3` and `pieces = "packed"`; legend
  position and margins through `ggplot2::theme()`. `title_gap` sets how far
  beyond its pebble a title box sits (default 0.06, as the pipeline's
  helper). Hyperedge labels
  (`edge_labels`) on a row of pieces sit inside their pebbles. New
  arguments: `titles` (`TRUE` for names, or a numeric selector for a count
  line), `title_prefix`, `notes` (a further line, from a named vector or a
  `group`/`note` table), `unit` (titles the colour legend and names counts,
  node area and transition widths), `pieces`. The caption and size legend
  call a node an "event".
  `node_sizes` without `direction` now draws points with an area legend
  (`scale_size_area()`); circles in data units are drawn with `direction`.
  `node_sizes`, `direction` and `transitions` read tables such as
  `(state, trials)` and `(from, to, trials)` as they are. The caption is one
  line ("Circle area: ... Triangle: points to the ... that most often
  follows it."). `tools` joins Imports (`toTitleCase()`).

* **`ring_sequences` and `ring_communities`**: 200 simulated walks with
  planted memory modules (four groups of four actions on a ring, shared
  actions between neighbours) and the planted community of every
  second-order state, built by `data-raw/ring_sequences.R`. They replace the
  simulator the `hon_communities()` walk-through defined inline.

* **`plot.net_hon_communities()` is rebuilt on it.** The default physical
  view is the community hypergraph (member = physical node, group =
  `"Community k"`) drawn by `plot.net_hypergraph()`: pebbles coloured by
  community (discrete Okabe-Ito, legend "Community"), in the blob look with a
  title box per community giving its name and flow,
  circles sized by physical flow, triangles pointing along the projected
  link flow, and a caption naming flow, one-node communities and the
  zero-flow states left out. It returns the ggplot; `...` reaches
  `plot.net_hypergraph()`. The custom circular layout, pie nodes and base
  legend are gone. `type = "states"` still draws with cograph.

* **`hon_outcome()` is withdrawn** (added in 0.5.0). Its per-actor features --
  the visit-weighted mean of a higher-order node's centrality or embedding
  over the actor's visits -- are a construction with no published basis, and
  nothing unreferenced ships in this package. The code is kept outside the
  package for later work. `sandwich` leaves Suggests with it.

* `hon_communities()` finds modules in a higher-order network with the map
  equation for memory networks (Rosvall et al. 2014; Edler, Bohlin & Rosvall
  2017). State nodes are clustered but coded over physical nodes, so a
  physical state can belong to several modules. Flow and codelength match the
  Python `infomap` package to ~1e-14 bits; the search (node aggregation with
  fine-tuning over seeded trials, trial stability as ARI) is reported beside a
  first-order map of the same flow, so the bits saved by memory are explicit.
  Returns `net_hon_communities`; `as.data.frame(what = "states" | "physical" |
  "modules" | "trials" | "first_order" | "codelength")`; `plot()` draws the
  physical view through `plot.net_hypergraph()` (see below) and uses cograph
  only for `type = "states"`. Coarse-tuning (Edler et al. 2017, Alg. 6) is not implemented.

* **`plot.net_markov_stability()` is redesigned** around four views chosen
  with `what =`. The default `"landscape"` places each state by stationary
  share (x, log scale) and persistence (y); dashed lines at the even share
  `1/n` and the mean persistence split the states into hubs, relays, traps
  and transients, point size is the mean stay, and memory states (`"a -> b"`)
  get their own shape. The other views delegate to cograph:
  `"states"` to `cograph::plot_centrality()`, `"passage_time"` to
  `cograph::plot_heatmap()` (states in share order, fastest and slowest
  passages named in the subtitle), and `"network"` to `cograph::plot_tna()`
  (node size = share, pie = persistence, the 10 most common states by
  default; returns `x` invisibly). `...` reaches the cograph function. Calls
  that pass `metrics` without `what` still get the per-metric view. To
  support the network view `markov_stability()` now also stores the
  row-normalised transition matrix as `$transition` (an additive field; no
  computed value changes). `grDevices`, already used by the hypergraph
  plots, is now declared in Imports.

* `as.data.frame.net_markov_stability()` gains `decreasing =` (so the
  shortest passage times and least persistent states are one call, not a
  subset) and `from =` / `to =` filters on the first-passage table. Ties
  break on the row key in both directions. Default output is unchanged.

* `markov_stability()` on a chain where pruning left a state with no outgoing
  transition now raises `hypernets_not_ergodic` (as well as
  `hypernets_bad_input`), and the message says to lower `min_freq` in
  `build_hon()`. The documentation shows the case on `ai_long`.

* Tests and examples call `build_hypa(order =)` instead of the deprecated
  `k =`, clearing the 45 deprecation warnings from the suite. The deprecated
  argument still works and still warns; a test now asserts both.

# hypernets 0.5.0

* `hon_outcome()` is the package's first verb that relates higher-order
  structure to an outcome. Per-actor features come from `hon_centrality()` and
  `build_honem()` -- no new mathematics -- aggregated as an exposure-weighted
  mean over the actor's own visits, decoded against the network by the
  longest-suffix rule. Returns a `net_outcome` whose `as.data.frame()` is one
  row per feature with `estimate`, `std_error`, `conf_low`, `conf_high`,
  `statistic`, `p` and `p_adj` (BH by default, and the result records which
  correction was applied). `nested_in =` switches to cluster-robust standard
  errors with t(G-1) intervals, and `vcov` selects the estimator: `"CR3"`
  (default, the cluster jackknife, `sandwich::vcovCL(type = "HC3")`) or
  `"CR1"` (Stata's, HC1). Measured interval coverage over 2000 replicates:
  0.9530 model-based, 0.9525 CR3, 0.9330 CR1 -- and 0.7200 when the
  clustering is ignored, which is why `nested_in` exists. CR3 is the default
  because CR1's shortfall is a known property of that estimator at moderate
  cluster counts, reproduced independently through `sandwich` itself.
  Collinearity, rank deficiency, few clusters, dropped actors and binomial
  separation all raise classed conditions; nothing is dropped silently.

* `markov_stability()` describes the random walk a transition matrix carries:
  persistence, stationary distribution, mean recurrence time, sojourn time and
  Kemeny-Snell mean first passage. On a `net_hon` the matrix is row-stochastic
  over *memory* states, so this is a higher-order random-walk analysis --
  stationary mass on memory states and first passage between them. A reducible
  chain now raises `hypernets_not_ergodic` instead of returning the 1e15
  artefacts an unguarded solve produces.

* **The memory verbs no longer read a long event table as a wide one.**
  `build_hon(long)` without `action`/`actor`/`time` used to treat each ROW as
  a trajectory, promoting actor ids and timestamps to states -- a two-actor,
  four-turn table became states `1, 2, 3, 4, A, B, C, s1, s2` and eight
  trajectories instead of three states and two, with no error and no warning.
  `build_hon()`, `build_mogen()`, `build_hypa()`, `markov_order_test()`,
  `bootstrap_hon()` and `compare_hon()` now raise `hypernets_long_format`
  (inheriting `hypernets_bad_input`) when handed a frame carrying two or more
  of the canonical long column names, with the remedy appropriate to that
  verb -- the long-format arguments where they exist, splitting the table
  where they do not.

* `hg_sequences()` closes the text-to-memory gap: a `text_hypergraph()` plus a
  `hg_cluster()` partition becomes the long `actor` / `time` / `action` table
  the memory family reads, so a transcript can go from topics to a
  higher-order transition model without the caller hand-writing a
  merge-order-split join. Note that `build_hon()` and friends need the
  `action` / `actor` / `time` arguments given explicitly.

* `as.data.frame()` gains a method for `net_markov_order_group`, which had a
  print method and no accessor. The accessor-contract test now DISCOVERS the
  package's `net_*` classes from `NAMESPACE` instead of comparing against a
  hand-maintained list, which is how that gap had stayed invisible.

* The local equivalence suite's gate now accepts the documented
  `HYPERNETS_EQUIV_TESTS`. It previously read only the pre-rename
  `HONETS_EQUIV_TESTS`, so the command in the project documentation skipped
  all 31 oracle files and reported "PASS 0" -- green, having run nothing.
  Cross-package identity against Nestimate's memory family is now proven from
  this side, in `local_testing_and_equivalence/test-identity-nestimate-memory.R`.

# hypernets 0.4.8

* `plot()` on a hypergraph lays each connected component out on its own and
  packs the results into a roughly square frame. A force-directed layout has
  no force between two disconnected components -- they repel and nothing pulls
  back -- so laying the whole hypergraph out at once let a stray hyperedge
  drift to the edge of the picture and set the scale for everything else. On a
  corpus of one connected core and three islands, the core was left 0.05 of
  the frame on average (0.02 at worst) over 20 seeds; it now gets 0.43 (0.24
  at worst). `padding` also sets the seam the packing leaves between
  components, so two components' pebbles never touch and imply a member they
  do not share.

* The force-directed layouts stop repelling beyond 2.5 ideal edge lengths --
  Fruchterman and Reingold's own grid variant. Unbounded repulsion left a
  weakly attached cluster no equilibrium: pushed by every vertex and pulled
  back by one edge, it settled roughly `n^(1/3)` ideal lengths out and the
  hyperedge bridging it stretched across the page. On a core-plus-satellite
  fixture over 30 seeds the satellite came in from 0.84 of the frame away to
  0.70, the core widened from 0.31 of the frame to 0.43, and the bridging
  pebble narrowed from 0.43 to 0.37. It also left fewer non-members inside a
  pebble than the uncut layout (0.20 per layout against 0.32).

* `plot(dismantled = TRUE)` no longer cuts a panel title to the panel's
  width, which silently turned `32006L0123` into `2006L012` -- an identifier
  that reads as a different hyperedge. Long names now overflow their strip;
  lower `edge_label_size` or `ncol` if they collide.

# hypernets 0.4.7

* `plot()` on a hypergraph draws every hyperedge as a smooth *pebble*: the
  convex hull of its members, widened by `padding` and low-pass filtered in
  the Fourier domain, then pushed out along its normal so every member keeps
  room inside. Overlapping pebbles are parted by white seams. New `detail`
  (smoothing; default `5`, `Inf` for the unsmoothed rounded hull) and
  `outline` (`"white"` seams, `"fill"` for an outline in the pebble's own
  colour, or any colour). Defaults are now `alpha = 0.45`, `linewidth = 1.1`
  and `padding = 0.045`. Drawing is plain ggplot2; `cograph::plot_simplicial()`
  is no longer used.

* The default `layout` is now `"bipartite"`: nodes and hyperedges placed
  together, so each hyperedge has a position of its own and its label sits
  there. The force-directed layouts (`"bipartite"`, `"spring"`) use a
  Fruchterman-Reingold placement computed in base R; on three test
  hypergraphs it left 0.3, 0.3 and 0 non-member nodes inside a pebble on
  average against 5, 1.8 and 10 for `cograph::layout_spring()`. `seed` no
  longer disturbs the caller's random number stream.

* New `center`: node names to place in the middle of the picture. Names
  absent from the hypergraph are skipped, so one vector serves several
  subsets. Centring costs some clarity, since convex pebbles around the
  hyperedges must cover the middle.

* `dismantled = TRUE` now returns one ggplot with a facet per hyperedge
  (members inked, others grey) instead of a gridExtra `gtable`, and is tested.

* Hyperedge labels have dark text in a box bordered with the hyperedge's
  colour; a label whose own position falls outside its pebble moves to the
  members' centroid. A whole-number colour scale (`color_by = "size"`) gets
  whole-number legend breaks. The same variable given to `color_by` and
  `linetype_by` yields one merged legend.

# hypernets 0.4.6

* The sparse-storage rule of `text_hypergraph()` (bag and sentence
  constructions) no longer overflows. It multiplied two integers
  (`n_docs * nrow(vocabulary)`), so any corpus large enough to need sparse
  storage overflowed `.Machine$integer.max` to `NA`, `isTRUE(NA)` was
  `FALSE`, and exactly the largest corpora went down the dense path and hit
  the vector memory limit inside `matrix()`. The product is now taken in
  double.

* `group_hypergraph()` refuses a dense incidence of more than
  `.Machine$integer.max` cells with the classed error
  `hypernets_dense_too_large`, instead of overflowing its flat cell index or
  attempting the allocation.

* `text_hypergraph()` gains `max_words` and `coverage`. With `min_count` they
  are one filter over one ranking (decreasing corpus count, ties broken
  alphabetically), so all three keep prefixes of the same order and the
  strictest wins. Both default to no pruning. Pruning emits a suppressible
  message and records `min_count`, `max_words`, `coverage`,
  `n_vocabulary_full` and `token_share` in the text layer. The `"knn"`
  construction rejects both, alongside the other token-based arguments.

# hypernets 0.4.5

* `plot()` on a hypergraph gains `dismantled` and `ncol`, drawing one panel
  per hyperedge on a shared layout, each titled with its hyperedge name.
  Overlaid blobs mislead: a blob is a hull around its members, so two
  hyperedges sharing nothing still overlap on the page wherever their hulls
  sweep past each other. Needs the suggested `gridExtra` and returns a
  `gtable`, so `edge_labels` does not apply. **This one has no unit tests
  yet.**

* `hg_cluster(what = "eigenvalues")` gains `n`, the number of leading rows to
  keep (default `Inf`, all of them), matching the argument `hg_centrality()`
  and `hg_pagerank()` already take. The spectrum carries one eigenvalue per
  node, so a corpus of a few thousand documents returned a few thousand rows
  when only the leading gaps decide how many groups the structure supports.

# hypernets 0.4.4

* `plot()` on a hypergraph gains `edge_labels` and `edge_label_size`, naming
  the hyperedges on the figure. A plot of the highest-ranked citation blocks
  was otherwise anonymous -- the ranking table names them and the picture did
  not -- and `labels`, which writes node names, is no help: a figure with
  seven hyperedges can carry three hundred nodes. Each label sits just outside
  the member furthest from the centre of the layout, clear of the overlap
  where the blobs meet, and the panel limits now expand to hold the labels so
  the outermost blob's label is not clipped away.

# hypernets 0.4.3

* New `hg_hypa()`: hypergeometric anomaly detection for co-occurring node
  pairs, the hypergraph counterpart of [build_hypa()]. A pair's propensity to
  share hyperedges is the product of the two hyperdegrees, and the observed
  co-occurrence count is referred to a hypergeometric law, so the test is
  analytic and needs no resampling. It returns one row per pair, where
  `hg_null_test()` returns one row for the whole hypergraph.

  Two things make it usable at a scale `build_hypa()` cannot reach. Only
  co-occurring pairs are scored, so the propensity stays sparse instead of
  the dense n-by-n outer product the memory family materialises; and
  `min_count` keeps pairs that cannot reach significance out of the
  multiplicity correction, chosen by count and never by p-value. On an EU
  citation hypergraph of 117,633 nodes and 181,364 hyperedges it scores
  435,784 pairs in 1.5 s, where `build_hypa()` on the same data exhausts
  32 GB.

  Results are ordered by adjusted significance rather than by `ratio`: ratio
  is maximised by the rarest pairs sitting just above `min_count`, which
  buries the heavily-cited pairs that actually depart from the null.

# hypernets 0.4.2

* `hg_edges()` no longer computes `n_neighbors` unless it is asked for.
  That measure takes two products over the node-by-node adjacency, which
  densifies on a hub-heavy network such as a citation graph; every other
  measure paid for it silently. On a temporal hypergraph, where `hg_edges()`
  loops over the event-time grid, the cost was paid once per snapshot: a
  directive co-citation hypergraph with 1,609 grid points now returns
  `measure = "n_incident_edges"` in 1.6 s. The parts of the edge table that
  do not depend on `s` are also computed once instead of once per value, so
  an `s` sweep no longer repeats them.

* `temporal_hypergraph()` and `group_hypergraph()` gain `separator`, which
  splits the `actor` column into one row per member before building.
  Bibliographic exports ship a hyperedge's members as a single delimited
  cell -- EUR-Lex `citationcelex` and `eurovoc`, Scopus and Web of Science
  reference and keyword fields -- and every caller was writing the same
  split, trim and drop-empties preamble. The incidence and edge tables are
  identical to the hand-rolled explode.

# hypernets 0.4.1

* `clean_text()` and `text_hypergraph()` gain `min_chars`, a minimum word
  length. Corpora extracted from PDFs carry single letters and short
  fragments -- the initials in "J.-P. Puissochet", enumeration markers, and
  the halves of words split by a hyphenated line break -- which survive stop
  lists and document-frequency floors because they are neither stop words nor
  rare. `text_hypergraph(min_chars = 3L)` gates the vocabulary itself, so no
  short token reaches a keyword table; `clean_text(min_chars = 3L)` does the
  same to the text and collapses the periods that removing initials leaves
  behind. Both default to keeping everything, so no existing result moves.
  `min_chars` joins `stop_words`, `min_count` and `weight` as an argument
  `construction = "knn"` refuses.

# hypernets 0.4.0

* The package is named **hypernets**, matching its repository and site.
  Every `honets_*` condition class is now `hypernets_*`, the package
  documentation is `?hypernets`, and `library(hypernets)` replaces
  `library(honets)`. No function, argument or result changed.

# hypernets 0.3.12

* `hg_subset()` gains `size =`, keeping hyperedges by member count, so a
  window snapshot can be cut to its three-code sets for `hg_motifs()`.

# hypernets 0.3.11

* `temporal_hypergraph()` reads a sequence table -- a wide data frame of
  states, a list of character vectors, or a `tna` / `netobject` model -- as
  one hyperedge per session whose states are contacts at their position, on
  a `"step"` clock: simple co-occurrence with time as order. In any input
  shape a membership may now carry its own time; the hyperedge spans its
  first to its last membership and a snapshot keeps the memberships present
  in its window, so `mode = "cumulative"` at step `t` is what a session had
  shown by `t` and `window = 3` is what it showed on three consecutive
  steps. `hg_growth()` counts the memberships present. Constant-time data
  are unchanged: every legal-hypergraphs number is identical.

# hypernets 0.3.10

* Two tests loaded `human_long` from the retired hypernets package instead of
  hypernets' own bundled copy (identical data), failing R CMD check on any
  machine without that package.

# hypernets 0.3.9

* The Legal hypergraphs workflow is a hand-knit document in `docs/`, not a
  vignette, so pkgdown no longer publishes it as an article. Its Figure 6d
  closeness now runs on the aggregate with duplicate tribunals collapsed,
  the paper's representation, and every chunk does one step.

# hypernets 0.3.8

* The pkgdown workflow installs igraph and gridExtra like the check workflow
  does, so the legal vignette's hyperedge betweenness (cograph -> igraph)
  renders on the site.

# hypernets 0.3.7

* hypernets installs without torch again. Five `torch::nn_module()` classes
  (HyperGAT layer and network, HGAT layer, AllSet DeepSets and PMA blocks)
  were built at namespace load, so `R CMD INSTALL` failed wherever the
  Suggests package torch was absent. They are now built on demand, and a
  regression test forbids any top-level reference to a Suggests package.

# hypernets 0.3.6

* pkgdown CI builds into `pkgdown/` (`dest_dir`), not `docs/`, which holds
  the hand-knit documents.

# hypernets 0.3.5

* The GitHub repository is `mohsaqr/hypernets`; `DESCRIPTION`, the README
  install line and the pkgdown site URL now point there. The package name is
  unchanged.

# hypernets 0.3.4

The temporal hypergraph speaks Dynet's vocabulary (`../temporal`), so a
relational log reads the same way in both packages.

* **`group`** replaces `cooccur_by` in `temporal_hypergraph()` and
  `group_hypergraph()` (Dynet's co-presence format is `actor`/`group`).
  `cooccur_by` and `member` remain as deprecated aliases for one release and
  warn with a `hypernets_deprecated` condition.
* **Column detection** uses Dynet's alias table, case-insensitively:
  `Sender`/`Receiver`, `source`/`target`, `onset`/`terminus`, `timestamp`
  are understood without being named; an explicit name must exist as
  written (`hypernets_missing_column`).
* **Time parsing** follows Dynet: numeric times stay as they are (unit
  `"step"`); `Date`, `POSIXct` and character date-times become elapsed
  time since the earliest time in a unit chosen for the span (`"days"`
  beyond three days) or given as `time_unit`. The object stores
  `time_unit` and `origin`, `print()` reports them, every `time` column of
  a result is on that clock, and dates passed to `at`, `start`, `end` or
  the observation bounds are converted (or refused on a numeric clock).
  Raw character dates are never compared.
* **`time` is a contact clock.** A hyperedge with a `time` is an
  instantaneous event, as in Dynet's contact format; the former "growing"
  reading is `mode = "cumulative"` in `hypergraph_snapshot()`,
  `hypergraph_snapshots()`, `hg_growth()` and `hg_edges()`, and `print()`
  says so. `mode = "all"` is deprecated: it is `"cumulative"` with `at`
  unset. `evolution` is replaced by `format` (`"interval"` or `"contact"`).
* **`observation_start` / `observation_end`** with Dynet's meaning: they
  bound the snapshot times and the measurement grid
  (`hypernets_outside_observation`), an open-ended hyperedge is active
  through the end of observation, and the stored memberships are never
  rewritten.
* **The measurement grid**: `hypergraph_snapshots(x, start, end, step,
  window, mode)` -- `step` is how often to look, `window` how much time
  each look covers (`0` a point, `"all"` the whole period), `at` names
  instants -- and `hg_growth()` and `hg_edges()` take the same four.
  Without `step` and `at` the grid is the event times, as before.
* `as.data.frame(x, what = "memberships" | "edges" | "nodes")` for a
  temporal hypergraph; series plots draw the calendar when the hypergraph
  has one.
* Behaviour-neutral on the bundled legal data: Table 2, the component
  sweep and the motif counts of the Legal Hypergraphs vignette are
  unchanged.

# hypernets 0.3.3

Every line of code a user reads now follows the one-call-per-line rule: no
verb nested inside another call, no `as.data.frame()` inside a call, no `$`
reach into a `net_*` result -- in the vignettes, the roxygen examples, the
hand-knit `docs/` documents and the tests alike. Three small additions made
the tidy versions possible:

* **`hg_agreement()`** gains an `aligned` column in its summary (the nodes
  that stay with the majority of their `x` label in `y`; the sum of
  `overlap` over `what = "mapping"`) and `node` / `label` column selectors,
  one name for both labelings or two for `x` and `y`, so a classifier can
  be scored against a column of the corpus table directly:
  `hg_agreement(predictions, corpus, node = c("node", "doc"), label =
  c("predicted", "year"), what = "table")`.
* **`group_hypergraph()`** keeps every column that is constant within a
  hyperedge as an attribute in `edge_data`, as `temporal_hypergraph()` does,
  so `hg_subset(where =)` and `edge_source =` work on it without assembling
  the table by hand. A plain member/group table keeps its original layout.
* **`as.data.frame(hg, what = "nodes")`** for any `net_hypergraph`: one row
  per node with `degree`, plus `block` for the stochastic block model
  generator; `sort_by = "degree"` and `top` apply.

# hypernets 0.3.2

The Legal Hypergraphs workflow (Coupette, Hartung & Katz 2024) is now
reproduced figure by figure on the authors' released data, and the verbs
that were missing for it are in the hypergraph family.

* **`plot()` for `net_hypergraph`**: nodes on a spring (or circle) layout of
  the clique projection, every hyperedge as a smooth translucent blob around
  its members, drawn by `cograph::plot_simplicial()`; hypernets maps
  `color_by` and `linetype_by` (`"size"`, a column of the edge metadata, or
  a value per hyperedge) to blob colours and line types and adds the legend;
  `labels` renames nodes; a `layout` table can be reused across panels.
* **`hg_subset()`**: sub-hypergraph by hyperedge names, by a node set (the
  induced sub-hypergraph) or by hyperedge `source` (the hypergraph of one
  citing decision), keeping edge metadata and multiplicities, sparse or dense.
* **`hg_growth()`**: node, hyperedge, distinct-set and membership counts at
  every event time of a temporal hypergraph, with cumulative columns for
  interval data and `components = TRUE` for the number of components, the
  share of the largest and its diameter; returns a `hypernets_series` table
  with a `plot()` method.
* **`hg_representations()`**: the paper's Table 2, the same data as binary
  graph, multi-graph, binary hypergraph and multi-hypergraph with edge counts
  and degree statistics, for the clique or the citation graph.
* **`hg_compare_communities()`**: several `hg_communities()` fits side by
  side, with summaries, pairwise AMI/ARI/NMI, cluster-size distributions
  and, given the hypergraph, each medoid's quality on its own projection;
  `plot()` draws Figure 8b and, through `cograph::plot_heatmap()`, the
  AMI/ARI matrix of Figure 8c (`as.data.frame(what = "matrix")`).
* **`temporal_hypergraph()` is defined like a network.** It takes an edge
  list (`from`, `to`) or co-occurrence data (`actor`, `cooccur_by`), a
  clock (`time`, or `start` and `end`; the evolution follows), a node
  universe (`nodes`: names, or a table whose first column is the node and
  whose `start` column its entry time) and `sparse = TRUE`. Every other
  column constant within a hyperedge is kept as a hyperedge attribute (a
  `source` column feeds self-association). The former wide `member = c(...)`,
  `edge`, `source`, `attributes` and `evolution` arguments are gone.
  `group_hypergraph()` takes the same `actor`, `cooccur_by`, `from` and `to`
  names (`member` and `group` still accepted). An empty snapshot keeps the
  universe; `summary()` reports the observation window.
* `hg_edges()` evaluates temporal hypergraphs snapshot by snapshot (`at`,
  `snapshot_mode`, `multiedges`), takes several `s` thresholds, and adds
  `what = "summary"`; distribution tables are `hypernets_distribution` objects
  whose `plot()` draws the CCDF, one curve per date or threshold.
* `hg_measures()` adds `n_neighbors` to the node table and
  `what = "distribution"` and `what = "components"`.
* `hg_project()` adds `method = "citation"` (source-to-member graph,
  directed or not); `hg_communities()` and `hg_community_quality()` take
  `method = "citation"`, and Infomap can run with directed flow.
* `hg_null_test()` adds the statistics `repeated_edges` and
  `repeated_pairs` and the paper's degree-ordered `method = "assignment"`.
* `hg_motifs()` returns a `hypernets_motifs` table that keeps its null draws
  (`as.data.frame(what = "draws")`) and plots the null distribution with the
  observed count; the motif census is vectorised (about seven times faster
  per null draw, identical counts).
* `group_hypergraph()` reads sparse incidence in one pass, so the GFCC
  aggregate builds in under a second instead of twenty.
* **Datasets** `icsid_tribunals` (2,226 tribunal seats of 742 ICSID cases),
  `gfcc_decisions` (3,618 decisions) and `gfcc_citations` (77,284 citations
  in 46,257 blocks), rebuilt from the authors' Zenodo archive by
  `data-raw/legal_hypergraphs.R` (CC BY-NC 4.0, attribution in the help
  pages). `hg_subset(where =)` selects hyperedges by any attribute, and
  `edge_source` may name an attribute column (`"citing"`).
* The `legal-hypergraphs` vignette is rewritten on those datasets: Tables 1
  and 2, Figures 3 to 8 and the repeated-collaboration test.

# hypernets 0.3.1

* **Topic summaries with plots**: `hg_topic_sizes()` (document and
  weighted shares), `hg_topic_quality()` (NPMI coherence and share
  exclusivity of each topic's top words) and `hg_membership()` (fuzzy
  c-means membership of every document in every topic, from the spectral
  embedding). Each returns a data.frame with a `plot()` method.
* **`hg_relations()`**: the topic-by-topic co-occurrence network through
  shared vocabulary (full counting; `similarity =` association, cosine,
  jaccard, inclusion or equivalence), as a `source, target, weight` edge
  list or a `cograph_network`.
* `hg_keywords()` takes an external score table through `scores =` and a
  sentence hypergraph as `hg` for sentence scope; the long form carries
  `size`; the print method is compact; centrality defaults to PageRank.
  `text_hypergraph()` chooses sparse storage automatically.
* `hg_cluster(edge_weights =)`: hyperedge weights for either Laplacian,
  numeric or `"idf"` (each word weighted by its inverse document
  frequency). The default `type = "zhou"` cut reads only hyperedge
  membership, so `weight = "tfidf"` alone does not change it; this argument
  is how tf-idf reaches the clustering. Comparison:
  `docs/covid-weighting.html`.
* **`clean_text()`**: corpus repair before `text_hypergraph()` -- HTML,
  mojibake, citations and numbering, URLs and DOIs, copyright notices,
  numbers, custom `remove` patterns, a content floor -- returning the same
  rows so nothing drops silently. Vignette: `vignette("covid-topics")`,
  sixteen topics of the COVID-19 education literature in five calls.
* **Sentence hyperedges.** `text_hypergraph(construction = "sentence")`
  binds the words of each sentence (HyperGAT's construction over a whole
  corpus), with `as.data.frame(hg, what = "sentences")`.
  `hg_keywords(type = "sentence_centrality", sentences = )` ranks a topic's
  words by centrality among the topic's sentences.
* **Topic descriptions.** `hg_keywords()` gains `type =`: `"mass"` (the
  previous score, default), `"frequency"` (raw counts), `"ctfidf"`
  (Grootendorst 2022 class-based tf-idf, matched to BERTopic's
  `ClassTfidfTransformer` to 1e-12), `"centrality"` (the word's
  `hypergraph_centrality()` in the cluster's own word hypergraph, measure
  chosen with `centrality =`) and `"attention"` (summed HyperGAT word
  attention). `hg_hypergat(what = "attention")` returns that per-document,
  per-word attention table. `type` takes several scores at once; the table
  (class `hypernets_keywords`, new leading `type` column) has a `plot()`
  method: one panel per topic and score, bars of the score per word
  (`value = "share"` to show shares). `sort_by = "share"` ranks a topic's
  words by the fraction of their total score it holds (distinctive rather
  than heavy vocabulary) and `min_docs` is the support floor; the table
  gains an `n_docs` column, and the collapsed table a `size` column.
  `hg_agreement(what = "mapping")` maps each cluster of one partition to
  the cluster of another that holds most of it. Worked examples: `docs/levebee-topics.html`
  and `docs/covid-topics.html`.
* **BERTopic benchmark.** `benchmarks/run_bertopic_benchmark.R` runs the
  actual `bertopic.BERTopic` package (ten seeds, three variants) against
  `hg_cluster()` on R8, on both the tf-idf hypergraph and a kNN hypergraph
  built from the same sentence embeddings; paired effects with bootstrap CIs
  are reported in the benchmarks article.
# hypernets 0.3.0

* **The complete Hayashi clustering family is implemented.**
  `hypergraph_cluster(algorithm = "symnmf")` adds Algorithm 2, RDC-Sym,
  alongside the existing RDC-Spec implementation. New
  `hypergraph_joint_cluster()` implements the patent experiment's J-NMF
  (Eq. 18) and JS-NMF (Eq. 19) objectives when an auxiliary node-relation
  matrix is available. All NMF paths expose their objective trace,
  convergence state, iteration count, chosen restart and fitted factors.
* **Full HyperGAT semantic hyperedges.** `hg_hypergat(semantic = "lda")`
  now implements the Ding et al. (2020) LDA path in native R: online
  variational Bayes on labeled training documents, class-count topics, and
  per-document topic edges from the top words. Precomputed keyword lists are
  accepted for exact replay of official preprocessing artifacts; the default
  remains the backward-compatible sentence-only ablation.
* **The neural paper set is complete.** New `hypergraph_hypergcn()` provides
  dynamic, fast and one-edge HyperGCN; `hypergraph_hnhn()` exposes the HNHN
  alpha/beta normalization exponents; and `hypergraph_allset()` implements
  both AllDeepSets and AllSetTransformer. `heterogeneous_hgat()` implements
  Linmei et al.'s distinct heterogeneous document/topic/entity graph model
  with node- and type-level attention. Each method has direct equation or
  construction invariants plus end-to-end torch tests.

* **cograph is the shared graph and plotting engine.** It is promoted from
  Suggests to Imports: hypernets owns hypergraph construction, incidence algebra
  and hypergraph-specific transformations, then uses cograph for ordinary
  graph algorithms and rendering. Dynet remains a separate peer and is not a
  dependency.
* **Dual hypergraph verb names.** Descriptive names such as
  `hypergraph_edges()`, `hypergraph_project()`, `hypergraph_pagerank()` and
  `hypergraph_classify()` are exported as direct bindings to their compact
  `hg_*()` forms; no implementation is duplicated. The raw-text attention
  model is also available as `text_hypergat()`. Existing engine names
  `hypergraph_measures()`, `hypergraph_centrality()` and
  `hypergraph_cluster()` keep their established meanings.
* **The full *Legal hypergraphs* method layer is implemented.** New temporal
  hypergraphs and snapshots; binary/multi and self-association projections;
  s-betweenness and s-closeness of hyperedges; log-subhypergraph centrality;
  induced Y/T/O motif censuses with the HypergraphX configuration MCMC;
  repeated Infomap with AMI-medoid selection; ARI/AMI/NMI agreement; and
  coverage, weighted coverage, performance and modularity. The authors'
  released ICSID data reproduce the published aggregate motif counts exactly:
  Y = 478, T = 7, O = 0.
* **The full Legal Hypergraphs archive reproduction now passes.** A standalone
  runner rebuilds the sparse 3,618-node GFCC citation-block hypergraph,
  verifies all four association projections and all eight archived Figure 8
  medoids, and reproduces the active/aggregate ICSID Figure 6 centrality
  rankings and values. `group_hypergraph()` now accepts an explicit node
  universe and sparse incidence, and normalized hyperedge centrality follows
  the exact NetworkX/HypergraphX convention on disconnected line graphs.

## The text family: texthypergraph folds into hypernets

The `texthypergraph` package is retired and its whole surface now lives here,
with its history of tests intact (491 shipped expectations came across; the
merged suite passes in full). hypernets gains a fourth family, **text
hypergraphs**, and the hypergraph family gains every generic method that
texthypergraph had built under its frozen-Nestimate contract.

* **New family — text hypergraphs.** `text_hypergraph()` (bag of words with
  smoothed tf-idf, token windows, or embedding kNN; dense or sparse),
  `stop_words_en()`, and the tidy verbs `hg_measures()`, `hg_centrality()`,
  `hg_cluster()`, `hg_keywords()`, `hg_classify()`, `hg_stability()`,
  `hg_agreement()`, `hg_seeds()`. Data: `covid_abstracts`,
  `covid_embeddings`. Vignettes `text-hypergraphs` and `text-constructions`.
* **Hypergraph family additions.** `knn_hypergraph()`, `dual_hypergraph()`,
  `hg_pagerank()` (personalized EDVW PageRank, sparse-capable),
  `hg_project()` (clique and association weightings), `hg_line_graph()`,
  `hg_edges()`, `hg_null_test()` (swap and configuration nulls), and the
  neural tier `hg_neural()` (HGNN) and `hg_hypergat()` (HyperGAT) in native
  torch. Sparse `Matrix::dgCMatrix` incidences are supported end to end by
  `hg_cluster()`, `hg_classify()`, `hg_pagerank()` and `hg_measures()`.
* **`hypergraph_transduction(normalization = )`.** New argument, default
  `"none"` (the raw Zhou 2006 argmax, unchanged behaviour). `"class_mass"`
  divides each class column by its total spread mass before the argmax (Zhu,
  Ghahramani & Lafferty 2003); without it, class-imbalanced seeds collapse
  every prediction onto the majority class (on R8, every test document is
  predicted "earn"). The result object records `$normalization`.
* **Nestimate dependency removed.** texthypergraph delegated its incidence
  construction and measures to Nestimate; those engines already lived here,
  and `group_hypergraph()` is `identical()` to `Nestimate::bipartite_groups()`
  up to the `member`/`player` argument name (tested). The package now
  imports only cograph, ggplot2, graphics, grid, Matrix, methods, parallel,
  RSpectra, stats and utils.
* **Collisions resolved without a value change.** texthypergraph carried a
  verbatim copy of the spectral trio; hypernets' copies were kept (they carry
  the plot methods, `top =`, and scalar `edge_weights`), and only the
  normalization argument was ported. `hg_pagerank()` agrees with
  `hypergraph_centrality(type = "pagerank")` to `1e-10` (tested);
  `hg_project(method = "clique")` equals `clique_expansion()` (tested);
  `text_hypergraph(construction = "window")` and `window_hypergraph()`
  produce the same incidence matrix on their shared domain (tested).
* **Condition classes** of the incoming code are `hypernets_*`
  (`hypernets_bad_input`, `hypernets_no_converge`,
  `hypernets_hypergraph_disconnected`, `hypernets_empty_corpus`,
  `hypernets_dropped_documents`, `hypernets_missing_embeddings`,
  `hypernets_missing_torch`, `hypernets_nonpositive_similarity`,
  `hypernets_sparse_unsupported`, `hypernets_sparse_too_large`,
  `hypernets_configuration_collapse`).
* **Infrastructure.** GitHub Actions (`R-CMD-check`, `pkgdown`), a pkgdown
  reference index covering every topic, `VignetteBuilder: knitr`, the
  benchmark harness under `benchmarks/` (build-ignored), and the
  text-family equivalence suites (HyperNetX, HyperG, DHG, the official
  HyperGAT code) under `local_testing_and_equivalence/`.

# hypernets 0.2.1

## Every accessor takes `top =`

Every verb that can return many rows now takes a `top =` argument, so a
caller never has to write `head()` around a result:

```r
as.data.frame(mo, what = "transitions", order = 2, top = 4)   # not head(..., 4)
```

`top` is applied **last** - after `what`, after every filter (`order_min`,
`min_count`, `dim`, `k`, `dimension`, `significant`), and after `sort_by` -
so `sort_by` and `top` compose: `top = n` is the first `n` rows of the table
as ordered. `top = NULL` (the default) returns everything, and the default
return of every accessor is unchanged. Semantics match the `top` that already
shipped on `path_counts()` and `pathways()`.

Added to: `as.data.frame()` for `net_hon`, `net_honem`, `net_hypa`,
`net_mogen`, `net_markov_order`, `net_path_dependence`, `net_hon_boot`,
`net_hon_compare`, `net_simplicial`, `net_q_analysis`,
`net_persistent_homology`, `net_persistence_landscape`, `net_hypergraph`,
`net_hypergraph_measures`, `net_hypergraph_cluster`,
`net_hypergraph_transduction`; and to `mogen_transitions()`,
`hon_centrality()`, `simplicial_degree()`, `hypergraph_centrality()`.

`path_counts(top =)` now validates its argument like the rest of the family:
a non-whole value such as `top = 2.5` is an error rather than a silent
truncation to 2, matching how `k` already behaved.

## Bug fix

* `plot.net_path_dependence()` clipped the label of its highest-KL context -
  the row the plot exists to show. The modal-flip labels are drawn to the
  right of each point, and the panel did not extend past the largest value,
  so ggplot cut the label off. The x scale now leaves room for it.

## Documentation

`docs/` (build-ignored) is reorganised from 23 per-verb vignettes into **four
documents**, one per structure family plus an overview:
`overview`, `memory-networks`, `simplicial-complexes`, `hypergraphs`. Each
section is one verb, worked end to end on the same data; the per-verb sources
are archived under `docs/_sections/`.

They gain **37 figures**, drawn with cograph. hypernets results are dual-classed
`cograph_network`, so `cograph::splot()` and `cograph::plot_simplicial()`
take them with no conversion step. `cograph` is added to `Suggests` and every
plot chunk is guarded on it.

The prose was also rewritten to remove 37 `head()` and 12 `subset()` calls
that subset a returned table on the public surface - the idiom `top =` now
replaces.

# hypernets 0.2.0

hypernets becomes **the** higher-order networks package: one package covering
all three structure families of the higher-order literature (Battiston et al.
2020) — memory networks, simplicial complexes, and hypergraphs — under a
single taxonomy.

## Absorbed families

* **Simplicial complexes and topological data analysis**, moved verbatim
  from Nestimate 0.9.0: `build_simplicial()` (clique, Vietoris-Rips and
  pathway complexes), `betti_numbers()`, `euler_characteristic()`,
  `persistent_homology()`, `persistence_landscape()`,
  `bottleneck_distance()`, `simplicial_degree()`, `q_analysis()`,
  `verify_simplicial()`.
* **Hypergraphs**, moved verbatim from Nestimate 0.9.0 by way of the
  short-lived earlier `hypernets` scaffold (0.1.2, never released), which is folded in
  and retired: `build_hypergraph()`, `window_hypergraph()`,
  `group_hypergraph()`, `hypergraph_measures()`, `hypergraph_centrality()`,
  `hypergraph_laplacian()`, `hypergraph_cluster()`,
  `hypergraph_transduction()`, `clique_expansion()`.
* The consolidation **deletes 223 lines of duplication** from hypernets'
  323-line `utils.R` -- 171 of them the clique-enumeration closure copied out
  of Nestimate's `simplicial.R`, the rest a `build_simplicial()` shim and a
  second copy of hypernets' own `.extract_edges_from_matrix()`. hypernets had to
  carry that closure because `build_hypergraph()` needs clique enumeration;
  with both families in one package, `build_hypergraph()` calls the real
  `build_simplicial()` again. The seven copied helpers were verified
  byte-identical to Nestimate's originals before the copy was removed.

## Taxonomy: renames

All renames are to the **surface only**. No computed value changed anywhere —
the `identical()` contracts against Nestimate 0.9.0 still hold for both
absorbed families, with only these names normalised away
(`local_testing_and_equivalence/test-identity-nestimate-{hypergraph,simplicial}.R`).

Every result class now carries a `net_*` class:

| Was | Is |
|---|---|
| `simplicial_complex` | `net_simplicial` |
| `persistent_homology` (class) | `net_persistent_homology` |
| `q_analysis` (class) | `net_q_analysis` |
| `persistence_landscape` (class) | `net_persistence_landscape` |
| `hypergraph_measures` (class) | `net_hypergraph_measures` |

Constructors and arguments:

| Was | Is | Why |
|---|---|---|
| `bipartite_groups()` | `group_hypergraph()` | it returns a hypergraph, not bipartite groups; now matches `window_hypergraph()` |
| `bipartite_groups(player =)` | `group_hypergraph(member =)` | the column names hypergraph nodes, which need not be people |
| `build_hypergraph(method =)` | `build_hypergraph(type =)` | same construction axis as `build_simplicial(type =)` |
| `params$method` | `params$type`, plus `params$source` | every hypergraph constructor now records `source`, so a result says how it was built |

Classed conditions are now uniformly `hypernets_*`:
`hypernets_no_converge` and `nestimate_hypergraph_disconnected` became
`hypernets_no_converge` (shared with the memory family's power iteration) and
`hypernets_hypergraph_disconnected`.

## Taxonomy: complete tidy-accessor coverage

Every `net_*` result class now has an `as.data.frame()` method, so no result
object requires reaching in with `$`. Eleven are new, each with a `what =`
argument for its secondary table:

* `net_hon` (`"rules"` / `"nodes"`, plus `order_min` and `sort_by`)
* `net_honem` (`"embeddings"` / `"variance"`)
* `net_hypa` (`"scores"` / `"over"` / `"under"`, plus `sort_by`)
* `net_mogen` (`"orders"` / `"transitions"`)
* `net_markov_order` (`"orders"` / `"null"`)
* `net_path_dependence` (plus `min_count`, `sort_by`)
* `net_simplicial` (`"simplices"` / `"f_vector"`, plus `dim`)
* `net_q_analysis` (`"q_levels"` / `"nodes"`)
* `net_persistent_homology` (`"persistence"` / `"betti"`, plus `dimension`,
  `sort_by`)
* `net_persistence_landscape` (plus `k`)
* `net_hypergraph_measures` (`"nodes"` / `"edges"` / `"global"`, plus
  `sort_by`)

A regression test asserts the coverage, so a future `net_*` class without an
accessor fails the suite.

## Structure

* `R/` is organised by family: `memory_*.R` (8 files), `simplicial_*.R` (4),
  `hypergraph_*.R` (7), plus shared `utils.R`, `pathways.R`, `data.R` and the
  package doc. Test files follow the same names.
* Nestimate's 1,561-line `simplicial.R` was split along its real dependency
  seams into `simplicial.R` (construction and structural measures),
  `simplicial_filtration.R` (the filtration and Z/2 boundary-reduction layer
  that both `build_simplicial(type = "vr")` and `persistent_homology()` sit
  on) and `simplicial_homology.R`. Code unchanged; the split was verified
  line-for-line content-preserving.

## Other

* `group_hypergraph()`'s weighted branch and `build_hypergraph()`'s incidence
  fill are vectorised (they were `for` loops accumulating into a matrix).
  Duplicate `(member, group)` cells are summed before assignment, which
  index assignment alone would not do.
* Package-level documentation (`?hypernets`) now states the three-family
  taxonomy, the verb grammar, and how the families cross into one another.

# hypernets 0.1.5

* New verb `hon_centrality()` (roadmap item A2): PageRank, betweenness
  and closeness computed on the higher-order topology and projected back
  onto first-order states (Scholtes, Wider & Garas 2016). Semantics
  follow pathpy 2.2.0 - that paper's reference implementation -
  generalized from fixed-order to the variable-order networks
  `build_hon()` produces, and verified against it: betweenness and
  closeness match exactly (< 1e-10) on second- and third-order
  topologies, PageRank to pathpy's own `tol = 1e-6`. The underlying
  kernels additionally match `igraph` on the same topologies.
  `project = FALSE` reports the centralities of the memory contexts
  themselves; `projection =` chooses how a higher-order node's PageRank
  is distributed ("scaled", "last", "first", "all"); `sort_by =` returns
  the table ranked.
* `build_hon()` gained the long-format interface already used by the
  inference verbs: `action`, `actor` and `time` column names, so an
  event table needs no manual splitting. Existing calls are unaffected
  (the arguments default to `NULL`) and produce byte-identical networks.
* New tutorial `Tutorial_docs/hon_centrality.html`: why the first-order
  network of a dense corpus cannot rank its states at all (complete
  digraph, uniform PageRank), what the higher-order ranking recovers,
  which contexts carry the flow, and when betweenness and closeness are
  saturated.

# hypernets 0.1.4

* New inference verbs for higher-order rules (roadmap item A1):
  `bootstrap_hon()` — sequence bootstrap with percentile CIs for rule
  probabilities and per-rule extraction *support*; `compare_hon()` —
  two-sample permutation comparison with per-edge BH adjustment and a
  pooled-count-weighted global test. Both precompute per-sequence counts
  once and rebuild replicates from reweighted counts (proven identical
  to re-counting the resampled multiset), draw all randomness serially
  (parallel runs reproduce serial results under a seed), accept long
  format (`action`/`actor`/`time`), and ship with `print`/`summary`/
  `plot`/`as.data.frame` methods (accessor filters `min_support`,
  `order_min`, `significant`, `sort_by`).
* Bundled example data `human_long` and `ai_long` (coded human-AI pair
  programming sessions, long format).
* New tutorial `Tutorial_docs/hon_inference.html`: rule stability by
  order, support-filtered reporting, and a diffuse early-vs-late shift
  (significant global test, no significant single edge).
* Internal: `.hon_extract_rules()` split into a counts-based core
  (`.hon_extract_rules_count()`); behavior unchanged (full equivalence
  suite re-verified).

# hypernets 0.1.3

* Roadmap: hypernets B2 (EDVW hypergraph PageRank) marked done in
  `EXPANSION-PLAN.md`. No package code changed.

# hypernets 0.1.2

* Roadmap: hypernets B1 (windowed sequence hyperedges) marked done in
  `EXPANSION-PLAN.md` with the shipped design recorded. No package code
  changed.

# hypernets 0.1.1

* Added the consolidated family expansion roadmap (`EXPANSION-PLAN.md`,
  build-ignored): hypernets higher-order features A1–A4 and the hypernets
  hypergraph sibling (scaffolded 2026-08-25). No package code changed.

# hypernets 0.1.0

* Initial release. Code moved from Nestimate 0.9.0 (delegation T0): `build_hon()`,
  `build_honem()`, `build_hypa()`, `build_mogen()`, `mogen_transitions()`,
  `path_counts()`, `markov_order_test()`, `path_dependence()`, and the
  `pathways()` generic with methods for `net_hon`, `net_hypa`, and `net_mogen`.
  Numbers are identical to the Nestimate implementations (same code, same RNG
  streams).
* Corrected the HONEM reference (Saebi, Ciampaglia, Kaplan & Chawla 2020,
  \doi{10.1089/big.2019.0169}); the author list previously cited was wrong.

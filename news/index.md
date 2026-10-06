# Changelog

## hypergraphs 0.6.9

- The Argentina example article gains a storyline of the busiest
  arbitrators and restores the centrality, per-community modularity and
  growth tables of its first version.
- The repository now holds only package files; development notes and
  local tooling are kept out of it.

## hypergraphs 0.6.8

- Test suite only: the regression test of default hull plots, whose
  baseline was recorded on macOS, is skipped on Linux and Windows, where
  the force-directed layout settles about 1e-4 away.

## hypergraphs 0.6.7

- [`hg_communities()`](https://mohsaqr.github.io/hypergraphs/reference/hg_communities.md)
  on a memory network numbers communities of equal flow the same way on
  every platform: flows equal to 12 significant digits are ties, broken
  by the first node.

## hypergraphs 0.6.6

- [`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)
  gains `collapse`: `FALSE` keeps every window as its own hyperedge, in
  order, named by its positions (`"1-3"`), with `sequence`, `start` and
  `end` in the edge metadata.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  hypergraph gains `type = "storyline"` for hyperedges with an order
  (the stored order or `sort_by`), such as the windows of one sequence:
  codes become lines that the windows gather.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a temporal
  hypergraph gains `type = "storyline"` (Tanahashi and Ma 2012): the
  `top` nodes with the most hyperedges (default 8) are lines from their
  first to their last hyperedge, hyperedges are columns in order of
  their start, and each column gathers the lines of its members. Lines
  are ordered by the barycentre rule (Sugiyama et al. 1981) over
  repeated sweeps, keeping the ordering with the fewest crossings.
  `start` and `end` limit the period; `edge_labels` and `point_size`
  style it. `spacing = "strength"` brings neighbouring lines closer the
  more hyperedges their nodes share (a full row for none, 0.35 for the
  strongest tie in the plot, linear in between) without changing their
  order. `width_by = "degree"` widens each line with its node’s number
  of hyperedges in the period, with a width legend. Each line has an
  Okabe-Ito colour and a point shape, distinct for up to 72 lines, named
  in the legend. `type` is now an explicit argument of
  [`plot.net_temporal_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/plot.net_temporal_hypergraph.md),
  and an argument of one view given to another raises
  `hypergraphs_bad_input`.

## hypergraphs 0.6.5

- Requires Nestimate (\>= 0.9.1), whose `prepare()` takes the `timezone`
  argument the sequence input passes; with 0.8.5 the vignettes failed.

## hypergraphs 0.6.4

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  hypergraph gains `type = "incidence"`, the incidence matrix after
  UpSet (Lex et al. 2014): one row per node ordered by hyperdegree, one
  column per hyperedge, a bar joining each hyperedge’s members,
  hyperdegree bars beside the rows and size bars above the columns when
  sizes differ. It stays legible where hulls overlap. The new `sort_by`
  orders the columns (`"size"`, an edge-metadata column such as a date,
  or one value per hyperedge; numbers largest first). `color_by`,
  `labels`, `edge_labels`, `label_size`, `edge_label_size` and
  `legend_title` apply; a hull-only argument raises
  `hypergraphs_bad_input`. Counts sort largest first; dates, text and
  clock columns (`start`, `end`, …) in time order.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  hypergraph sizes the node labels to the number of labelled nodes when
  `label_size` is not given: 4.2 mm up to 12 labels, then shrinking with
  the square root of the count to 2.2 mm from 44 labels on (was a fixed
  4.2 mm). An explicit `label_size` is used as given; an invalid one
  raises `hypergraphs_bad_input`.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a temporal
  hypergraph plots its snapshot as a hypergraph with
  [`plot.net_hg()`](https://mohsaqr.github.io/hypergraphs/reference/plot.net_hg.md)
  (hulls, or `type = "incidence"`), so the hyperedges are kept; it used
  to plot a pairwise projection through
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
  `method` is deprecated (`hypergraphs_deprecated`);
  `plot(pairwise_network(hg_snapshot(x, at)))` plots the projection.
- [`hg_snapshot()`](https://mohsaqr.github.io/hypergraphs/reference/hg_snapshot.md)
  keeps the data’s names for nodes and hyperedges, so a printed or
  plotted snapshot says “cases” and “arbitrator” instead of “edge” and
  “node”.
- [`hg_edge_centrality()`](https://mohsaqr.github.io/hypergraphs/reference/hg_edge_centrality.md)
  and
  [`hg_motifs()`](https://mohsaqr.github.io/hypergraphs/reference/hg_motifs.md)
  on a temporal hypergraph report the snapshot `time` as a date for a
  calendar hypergraph and as the number on the clock otherwise; it used
  to be a character label of the day offset (`"3104"`), also when `at`
  was given as dates.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  hypergraph with repeated member sets accepts `node_groups`: grouping
  the nodes now draws every hyperedge, as `color_by` does, instead of
  failing on the node sizes of the distinct-set view.
- `icsid_tribunals`: `respondent` and `subject` are trimmed of the
  spaces the source archive pads some values with. “Argentine Republic”
  held 12 of Argentina’s 47 cases, so a filter on the clean name missed
  them; there are now 145 distinct respondents instead of 201. No other
  column changes.
- `hg_subset(where =)` raises `hypergraphs_bad_input` for a value that
  no hyperedge takes, naming the closest values, instead of returning an
  empty hypergraph.

## hypergraphs 0.6.3

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  hypergraph gains `node_groups`: each node is coloured and shaped by
  its group, from a community fit
  ([`hg_communities()`](https://mohsaqr.github.io/hypergraphs/reference/hg_communities.md),
  [`hg_mmsbm()`](https://mohsaqr.github.io/hypergraphs/reference/hg_mmsbm.md),
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md))
  or any node/label table, so communities can be plotted over the
  hyperedges. The hulls turn grey unless `color_by` is given.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
  `color_by = "size"` lists only the sizes of the hyperedges it draws,
  and the curves of a distribution over calendar times are named by
  date.

- A clock column of a hyperedge (`start`, `end`, `time`, `session`, …)
  is no longer chosen as the count that titles and colours the hulls;
  the decision figures of a temporal snapshot are readable again.

- [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  on a hypergraph has a `weight` column only for a hypergraph of
  windows; other hypergraphs no longer print a column of `NA`.

- An unknown `what` in any accessor or verb raises
  `hypergraphs_bad_input` naming the available tables.

- [`hg_agreement()`](https://mohsaqr.github.io/hypergraphs/reference/hg_agreement.md)
  reads
  [`hg_mmsbm()`](https://mohsaqr.github.io/hypergraphs/reference/hg_mmsbm.md)
  and
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md)
  fits and text hypergraphs directly, and a table keyed by `doc` needs
  no `node =`.

- [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  on a community fit gains `converged`, `sort_by` and `top`.

- New article “Legal hypergraphs” on the tribunals of ICSID and the
  citations of the German Federal Constitutional Court, with hypergraph
  modularity.

- New dataset `forum_posts`: 360 simulated forum posts by 30 students
  whose sequence of topics has second-order memory, generated by
  `data-raw/forum_posts.R`. It is the example corpus for
  [`hg_sequences()`](https://mohsaqr.github.io/hypergraphs/reference/hg_sequences.md).

- `read_hif()` and `write_hif()` are removed, and with them the
  `jsonlite` suggestion. A hypergraph is exchanged with other packages
  as a `netobject` or a `cograph_network`.
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  on a hypergraph no longer offers `what = "node_data"` or
  `"incidence_data"`, which only a read HIF file filled.

## hypergraphs 0.6.2

- [`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md)
  gains `separator`: a delimited field such as the author keywords of a
  bibliographic export is read as whole-phrase terms, each a hyperedge
  over the papers that carry it, the keyword incidence of co-word
  analysis (Callon et al. 1983). Terms are lowercased and trimmed of the
  punctuation exports leave at their edges.

- New
  [`hg_dictionary()`](https://mohsaqr.github.io/hypergraphs/reference/hg_dictionary.md)
  labels each node by the dictionary category whose terms it carries
  most (Grimmer & Stewart 2013); the result is the `labels` input of
  [`hg_classify()`](https://mohsaqr.github.io/hypergraphs/reference/hg_classify.md)
  and
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md).

- Documents without a label are `split = "unlabelled"` in an
  `hg_classification`, so `hg_get(fit, split = "unlabelled")` reads the
  classifier’s proposals.

- [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  on a text hypergraph gains `node` (documents by id or by a table with
  a `node` column), `sort_by` and `top` (vocabulary).

- [`hg_classify()`](https://mohsaqr.github.io/hypergraphs/reference/hg_classify.md)
  and
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
  gain `holdout`: a share of the labels, drawn within each class with
  `seed`, is hidden, predicted and scored. The result is an
  `hg_classification` that prints the held-out accuracy and balanced
  accuracy (Brodersen et al. 2010).
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  reads its `"predictions"`, `"accuracy"`, `"classes"` and `"confusion"`
  tables, and filters the predictions with `split`, `correct`, `node`,
  `sort_by` and `top` (per class).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  confusion table. A
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
  result retains its trained network.
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  computes attention on request, and `plot(type = "attention")` shows
  its diagnostic weights. `plot(type = "hyperedges")` plots a document’s
  word hypergraph.

- `labels` of
  [`hg_classify()`](https://mohsaqr.github.io/hypergraphs/reference/hg_classify.md)
  and
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md),
  `clusters` of
  [`hg_keywords()`](https://mohsaqr.github.io/hypergraphs/reference/hg_keywords.md),
  [`topic_network()`](https://mohsaqr.github.io/hypergraphs/reference/topic_network.md)
  and
  [`hg_topic_sizes()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topic_sizes.md),
  and `group` of `hg_get(what = "prevalence")` accept the name of a
  document column (`labels = "subject"`).

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on an
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md)
  model gains `type = "prevalence"` with `group`, the prevalence of
  every topic in each group.

- [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
  now returns a reusable fitted classifier in every run.
  [`predict()`](https://rdrr.io/r/stats/predict.html) classifies new
  documents using its trained network, frozen vocabulary and topic
  keywords; optional observed labels evaluate these predictions without
  updating the model.
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  reads predictions, evaluation, document text, training history and
  optional attention diagnostics.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows
  confusion or training loss. The vocabulary is estimated from known
  labels only. Unscorable held-out and new labelled documents count as
  errors, and prediction preserves all new input rows.

- HyperGAT attention diagnostics include uniform normalization baselines
  and the fraction of words exclusive to each edge. Attention plots show
  the baseline alongside the learned weights. Attention is computed only
  when requested and is described as an internal weight, not a
  prediction explanation. Legacy `what` extraction remains available
  with a warning.

- [`hg_subset()`](https://mohsaqr.github.io/hypergraphs/reference/hg_subset.md)
  gains `component = "largest"`, which keeps the largest connected
  component, the input that label spreading and spectral clustering
  need. A subset text hypergraph now carries a matching text layer, so
  [`hg_keywords()`](https://mohsaqr.github.io/hypergraphs/reference/hg_keywords.md)
  and the other text verbs work on it.

- [`stop_words_en()`](https://mohsaqr.github.io/hypergraphs/reference/stop_words_en.md)
  gains `type = "snowball"`, the 174-word Snowball English stop list
  with pronouns and contractions, for chat and other informal text.

- [`clean_text()`](https://mohsaqr.github.io/hypergraphs/reference/clean_text.md)
  repairs emoji and other characters garbled by a UTF-8-as-Windows-1252
  export.

- Labels or groups that name documents
  [`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md)
  or
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
  dropped as empty are set aside with a warning in
  [`hg_keywords()`](https://mohsaqr.github.io/hypergraphs/reference/hg_keywords.md),
  [`hg_topic_sizes()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topic_sizes.md),
  [`topic_network()`](https://mohsaqr.github.io/hypergraphs/reference/topic_network.md),
  [`hg_classify()`](https://mohsaqr.github.io/hypergraphs/reference/hg_classify.md),
  [`hg_neural()`](https://mohsaqr.github.io/hypergraphs/reference/hg_neural.md)
  and
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md),
  so the corpus table can be passed back whole; an id that never was a
  document is still an error.

- [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  on an
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md)
  model gains `what = "prevalence"` with `group =`: the mean share of
  each topic within each group of documents, the topic prevalence by
  covariate of Roberts et al. (2014).

- [`clean_text()`](https://mohsaqr.github.io/hypergraphs/reference/clean_text.md)
  gains `boilerplate = TRUE`, which removes publisher names, company
  suffixes and the phrases of licence and rights notices wherever they
  occur in a text, so a classifier of bibliographic abstracts does not
  learn the publisher and the year from the notice.

- [`clean_text()`](https://mohsaqr.github.io/hypergraphs/reference/clean_text.md)
  removes the word “Copyright” together with the notice that follows it
  (“Copyright © 2020 Elsevier Ltd.” left “Copyright” behind).

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  [`simplicial()`](https://mohsaqr.github.io/hypergraphs/reference/simplicial.md)
  result gains `type =`. `"simplices"` (the default) draws the maximal
  simplices as before; `"summary"` draws the face counts, the Betti
  numbers and the simplicial degree, the summary figure of the complex.

- **Input formats with their own vocabulary.** Membership data name a
  `node` and a `hyperedge` (or `from` and `to`); event data name an
  `action` with its `session` and `actor`, as in the memory family;
  `group` is always the comparison variable, as in Nestimate. In
  [`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md)
  and
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/temporal_hypergraph.md),
  `group =` (the hyperedge column) is now `hyperedge =` and `by =` is
  now `group =`.
  [`hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/hypergraph.md)
  recognises the format from the arguments: `action` with `session` (or
  `actor`) gives one hyperedge per session, `window` gives windows,
  `node` with `hyperedge` gives membership hyperedges.

- **[`pairwise_network()`](https://mohsaqr.github.io/hypergraphs/reference/pairwise_network.md)
  builds the pairwise network of a hypergraph.** It returns a network
  object (`net_hg_pairwise`, also `netobject` and `cograph_network`) for
  every projection, chosen with `type = "clique"` (default),
  `"association"` or `"citation"`;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) plots it
  through cograph,
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  returns its edges (`what = "nodes"` its nodes), and
  [`hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/hypergraph.md)
  reads it as a network. It replaces `hg_project()` and
  `hg_clique_expansion()`, which are removed.

- **The topic network is
  [`topic_network()`](https://mohsaqr.github.io/hypergraphs/reference/topic_network.md)**,
  a constructor like
  [`pairwise_network()`](https://mohsaqr.github.io/hypergraphs/reference/pairwise_network.md).
  It replaces `hg_network()`, its name in 0.6.1, and `hg_relations()`,
  which are removed: `hg_relations(hg, clusters)` is
  `hg_topic_network(hg, clusters = clusters)`.

- **Constructors are bare nouns.**
  `random_hypergraph(type = "uniform" | "regular" | "gnp" | "sbm")`
  replaces `hg_sample_uniform()`, `hg_sample_regular()`,
  `hg_sample_gnp()` and `hg_sample_sbm()`, with the same models and the
  same seeded draws; its arguments are given by name. `read_hif()` and
  `write_hif()` replace `hg_read_hif()` and `hg_write_hif()`.

- **`group_hypergraph(min_share =)`** keeps the frequent sets that reach
  a minimum support, the share of a group’s sets that contain them
  (Agrawal & Srikant 1994); without `top`, every such set is kept. It
  applies to data frames, clustered sequences and topic models, and
  [`print()`](https://rdrr.io/r/base/print.html) states the rule.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  hypergraph packs disconnected pieces into one frame by default
  (`pieces = "packed"`), so each piece takes room in proportion to its
  size; `pieces = "row"` gives the former side-by-side frames.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  hypergraph whose hyperedges repeat (the trials of an event log) draws
  its distinct sets, each coloured by its number of copies; a hypergraph
  of counted sets with one group is drawn as that group; a window
  hypergraph is coloured by its window counts.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  [`pairwise_network()`](https://mohsaqr.github.io/hypergraphs/reference/pairwise_network.md)
  uses a circle layout. None of these needs an argument.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  hypergraph takes a numeric hyperedge attribute as its count only when
  the attribute varies between hyperedges; a constant one, such as the
  session number of the runs of one session, no longer colours and
  titles the hyperedges. `color_by = "weight"` colours the hyperedges of
  a
  [`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)
  by their window counts.

- **`tutoring_events` replaces `debug_events`.** It holds every step of
  the tutoring data, 84,356 events from 13,309 problem steps (20,626
  trials), with each event renamed to a word of the same meaning and
  similar events merged (19 events become 15). `debug_events` held a
  random 4,000 steps under a coding-assistant vocabulary. Its `outcome`
  column marks a step `completed` (with a correct answer) or `stopped`
  (without one).

- [`hg_measures()`](https://mohsaqr.github.io/hypergraphs/reference/hg_measures.md)
  builds the hyperedge-by-hyperedge overlap matrices only for
  `what = "overlap"`. The node table and the summary of a hypergraph
  with 20,000 hyperedges took minutes or exhausted memory; they take
  under a second.

- [`hg_subset()`](https://mohsaqr.github.io/hypergraphs/reference/hg_subset.md)
  keeps the window counts of a
  [`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)
  aligned with the hyperedges it keeps; printing or plotting such a
  subset failed.

- [`hg_agreement()`](https://mohsaqr.github.io/hypergraphs/reference/hg_agreement.md)
  accepts a fit of
  [`hg_communities()`](https://mohsaqr.github.io/hypergraphs/reference/hg_communities.md)
  and compares its medoid partition.

- New vignette, *Hypergraphs*: observed groups, frequent sets, windows,
  projection, node measures, centrality, communities and null models on
  `debug_events`.

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

- `hg_network()` replaces `hg_relations()` (kept as a deprecated alias
  that warns with `hypergraphs_deprecated` and returns the identical
  result). Besides a partition (`clusters =`), it takes a topic model
  (`topics =`). Without `threshold` two topics are linked by the
  correlation of their shares (stm’s simple `topicCorr()`, equal on the
  local oracle); with `threshold` the weight is the number of documents
  in which both topics reach that share (Abuhay et al. 2017; Cassi et
  al. 2017), normalisable by the same similarity measures.

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

#' hypergraphs: Higher-Order Network Analysis
#'
#' A higher-order network is a network in which a relation binds more than
#' two nodes at once or depends on more than the current node. The
#' literature describes three structures for such relations: memory
#' networks, simplicial complexes and hypergraphs. The package covers all
#' three, together with text hypergraphs, which represent a corpus as a
#' hypergraph. The works that define the structures are listed under
#' References.
#'
#' @section Families:
#'
#' \describe{
#'   \item{Memory networks}{A node is a state paired with the states visited
#'     just before it, so the probability of the next step depends on the
#'     path by which the current state was reached. The networks are
#'     estimated from categorical sequences given in long, wide, list or
#'     model form (see [sequence-input]). Constructors [hon()], [honem()],
#'     [mogen()] and [hypa()]; order selection [markov_order()], [memory()]
#'     and [hg_markov_stability()]; inference [hg_bootstrap()] and
#'     [hg_compare()]; measures [hg_centrality()] and [hg_communities()].}
#'   \item{Simplicial complexes}{A simplex is a set of nodes that are all
#'     related to each other, and a complex contains every subset of each of
#'     its simplices. Constructor [simplicial()]; measures [hg_betti()],
#'     [hg_euler()], [hg_degree()] and [hg_qanalysis()]; persistent homology
#'     [hg_homology()], [hg_landscape()], [hg_bottleneck()] and
#'     [hg_wasserstein()].}
#'   \item{Hypergraphs}{A hyperedge is an arbitrary set of nodes, and its
#'     subsets need not be hyperedges. The main constructor [hypergraph()]
#'     selects one of [group_hypergraph()], [window_hypergraph()],
#'     [network_hypergraph()] and [temporal_hypergraph()] from its input;
#'     [knn_hypergraph()], [dual_hypergraph()], [random_hypergraph()] and
#'     [hg_subset()] build further hypergraphs. Descriptive measures
#'     [hg_measures()], [hg_representations()], [hg_edges()],
#'     [hg_line_graph()], [hg_transitivity()], [hg_assortativity()],
#'     [hg_degree_correlation()] and [hg_motifs()]; centrality
#'     [hg_centrality()], [hg_pagerank()] and [hg_edge_centrality()];
#'     spectral methods [hg_laplacian()], [hg_cluster()] and [hg_classify()];
#'     communities [hg_communities()], [hg_modularity()], [hg_mmsbm()],
#'     [hg_community_quality()] and [hg_compare_communities()]; null models
#'     [hg_null_test()] and [hypa()]; time [hg_snapshot()], [hg_snapshots()]
#'     and [hg_growth()]; projection [pairwise_network()].}
#'   \item{Text hypergraphs}{A corpus is a bipartite structure of documents
#'     and words, and so a hypergraph in either orientation: documents as
#'     nodes joined by shared words, or words as nodes joined by shared
#'     documents. Constructor [text_hypergraph()] (bag of words, token
#'     windows, sentences or embedding nearest neighbours), with
#'     [clean_text()] and [stop_words_en()] for preparation; keywords
#'     [hg_keywords()]; topic models [hg_topics()] and [hg_topic_search()];
#'     clusters and their diagnostics [hg_cluster()], [hg_cocluster()],
#'     [hg_stability()], [hg_agreement()], [hg_membership()],
#'     [hg_topic_quality()] and [hg_topic_sizes()]; classification
#'     [hg_classify()] with [hg_seeds()]; topic relations [topic_network()];
#'     sequences of topics [hg_sequences()].}
#' }
#'
#' @section Naming:
#'
#' Constructors are nouns: [hon()], [simplicial()], [hypergraph()],
#' [text_hypergraph()] and the `<source>_hypergraph()` functions. Where a
#' family has several construction routes, `type =` selects one. All other
#' functions begin with `hg_`. A measure returns a data frame with one row
#' per node or per structure, and an inference function returns a result
#' object with estimates, intervals and p-values. [hg_centrality()],
#' [hg_communities()] and [hypa()] take a memory network or a hypergraph.
#' The tables of every result object are read with
#' `hg_get(x, what = )`, where `what` names a secondary table and the
#' named arguments `sort_by` and `top` order and shorten it.
#'
#' @section Links between families:
#'
#' A text hypergraph is a hypergraph, so every hypergraph function accepts
#' it. [simplicial()] with `type = "pathway"` builds a simplicial complex
#' from a memory network, and [window_hypergraph()] builds a hypergraph from
#' the same sequences. [pairwise_network()] projects a hypergraph to a
#' pairwise network. [hg_sequences()] turns a clustered corpus into
#' sequences for the memory family.
#'
#' @references
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
#' Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
#' \doi{10.1016/j.physrep.2020.05.004}
#'
#' Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University Press.
#' \doi{10.1017/9781108770996}
#'
#' @keywords internal
#' @importFrom cograph splot
#' @importFrom Nestimate build_simplicial
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom ggplot2 .data
"_PACKAGE"

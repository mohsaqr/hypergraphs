# hypergraphs

Higher-order network analysis in R.

A higher-order network is a network in which a relation binds more than two
nodes at once or depends on more than the current node. The literature
describes three structures for such relations, and hypergraphs covers all
three, together with text hypergraphs, which represent a corpus as a
hypergraph.

| Family | A relation is | Built from |
|---|---|---|
| Memory networks | a state paired with the states visited just before it | categorical sequences |
| Simplicial complexes | a set of nodes that are all related to each other, with every subset | a network, a distance matrix or a memory network |
| Hypergraphs | an arbitrary set of nodes | group memberships, sliding windows over sequences, network cliques, interval or contact records, or embedding neighbours |
| Text hypergraphs | a document and the words it contains, or a word and its documents | a corpus: bag of words, token windows, sentences or sentence embeddings |

The package depends on base R, ggplot2 and Matrix, with RSpectra for the
sparse spectral computations. Sentence embeddings need the optional `sbert`
package.

## Functions

### Memory networks

| Function | Method | Reference |
|---|---|---|
| `hon()` | Higher-order network with rule extraction (BuildHON+) | Xu, Wickramarathne & Chawla (2016) |
| `honem()` | Higher-order network embedding | Saebi, Ciampaglia, Kaplan & Chawla (2020) |
| `hypa()` | Hypergeometric path anomaly detection | LaRock et al. (2020) |
| `mogen()` | Multi-order generative model | Scholtes (2017) |
| `markov_order()` | Permutation test of the Markov order | Anderson & Goodman (1957) |
| `memory()` | Per-context divergence of order-k from order-1 transitions | Cover & Thomas (2006) |
| `hg_markov_stability()` | Persistence, stationary flow and first-passage times per state | Kemeny & Snell (1976) |
| `hg_bootstrap()` | Bootstrap confidence intervals and support of rules | Efron & Tibshirani (1993) |
| `hg_compare()` | Two-sample permutation comparison of rules | Good (2005) |
| `hg_centrality()` | PageRank, betweenness and closeness on the higher-order topology, projected to states | Scholtes, Wider & Garas (2016) |
| `hg_communities()` | Map equation for memory networks | Rosvall et al. (2014) |

### Simplicial complexes

| Function | Method | Reference |
|---|---|---|
| `simplicial()` | Clique, Vietoris-Rips or pathway complex, checked against the cliques of the graph and the Euler-Poincaré formula | Hatcher (2002) |
| `hg_betti()`, `hg_euler()` | Betti numbers and the Euler characteristic | Hatcher (2002) |
| `hg_homology()` | Persistence diagrams and Betti curves over a filtration | Edelsbrunner & Harer (2010) |
| `hg_landscape()`, `hg_bottleneck()`, `hg_wasserstein()` | Persistence landscapes and distances between persistence diagrams | Bubenik (2015); Cohen-Steiner et al. (2007); Kerber et al. (2017) |
| `hg_degree()`, `hg_qanalysis()` | Simplicial degree and q-connectivity | Atkin (1974) |

### Hypergraphs

| Function | Method | Reference |
|---|---|---|
| `hypergraph()` | Main constructor: selects one of the constructors below from its input | Battiston et al. (2020) |
| `group_hypergraph()` | Hyperedges from group memberships or from the most frequent sets | Perc et al. (2013) |
| `window_hypergraph()` | Hyperedges from sliding windows over sequences | Ding et al. (2020) |
| `network_hypergraph()` | Hyperedges from the cliques of a network | Burgio et al. (2020) |
| `temporal_hypergraph()`, `hg_snapshot()`, `hg_snapshots()`, `hg_growth()` | Interval and contact hypergraphs (`node`/`hyperedge` or `from`/`to`, calendar clocks, observation bounds), their active and cumulative snapshots, and growth over time | Coupette, Hartung & Katz (2024) |
| `knn_hypergraph()` | Hyperedges from the nearest neighbours in an embedding | Feng et al. (2019) |
| `dual_hypergraph()` | Exchange of the roles of nodes and hyperedges | Berge (1989) |
| `random_hypergraph()` | Uniform, regular, G(n, p) and stochastic block random hypergraphs | Marchette (2021) |
| `hg_subset()` | Sub-hypergraph by hyperedges, nodes, hyperedge attributes, size or component | Coupette, Hartung & Katz (2024) |
| `hg_measures()`, `hg_edges()` | Node, hyperedge and overlap measures, and hyperedge distributions | Lee et al. (2025); Do et al. (2020); Coupette, Hartung & Katz (2024) |
| `hg_representations()` | Graph and hypergraph representations of the same data, compared by size | Coupette, Hartung & Katz (2024) |
| `hg_transitivity()`, `hg_assortativity()`, `hg_degree_correlation()` | Clustering coefficients, degree assortativity and degree correlation | Klimm, Deane & Reinert (2021); Chodrow (2020); Lotito et al. (2023) |
| `hg_motifs()` | Census of Y, T and O motifs of three hyperedges, against a configuration model | Coupette, Hartung & Katz (2024) |
| `hg_centrality()`, `hg_pagerank()` | Clique, Z and H eigenvector centralities; PageRank of the random walk with edge-dependent vertex weights | Benson (2019); Chitra & Raphael (2019) |
| `hg_line_graph()`, `hg_edge_centrality()` | s-line graph and s-betweenness and s-closeness of hyperedges | Aksoy et al. (2020) |
| `pairwise_network()` | Pairwise projection with clique, association or citation weights | Coupette, Hartung & Katz (2024) |
| `hg_laplacian()`, `hg_cluster()`, `hg_classify()` | Hypergraph Laplacians, spectral and symmetric-NMF clustering, label spreading | Zhou, Huang & Schölkopf (2006); Hayashi et al. (2020) |
| `hg_communities()`, `hg_community_quality()`, `hg_compare_communities()` | Repeated Infomap with medoid selection, iteratively reweighted modularity maximization, and the quality and agreement of partitions | Coupette, Hartung & Katz (2024); Kumar et al. (2020) |
| `hg_modularity()` | Linear hypergraph modularity | Kamiński et al. (2019) |
| `hg_mmsbm()` | Mixed-membership stochastic block model | Ruggeri et al. (2023) |
| `hg_null_test()` | Swap, configuration and assignment null models | Gotelli (2000); Chodrow (2020); Coupette, Hartung & Katz (2024) |
| `hypa()` | Hypergeometric anomaly detection for co-occurring node pairs | LaRock et al. (2020) |

### Text hypergraphs

| Function | Method | Reference |
|---|---|---|
| `text_hypergraph()` | Corpus to weighted document-word hypergraph: bag of words with smoothed tf-idf, token windows, sentences or embedding nearest neighbours | Hayashi et al. (2020); Ding et al. (2020); Manning, Raghavan & Schütze (2008) |
| `clean_text()`, `stop_words_en()` | Text normalization and an English stop-word list | Manning, Raghavan & Schütze (2008) |
| `hg_cluster()`, `hg_keywords()` | Spectral or symmetric-NMF clustering of documents and the keywords of each cluster | Zhou, Huang & Schölkopf (2006); Hayashi et al. (2020) |
| `hg_topics()`, `hg_topic_search()` | Mixed-membership topic model (KL-divergence NMF) and the choice of the number of topics | Lee & Seung (2001); Roberts, Stewart & Tingley (2019) |
| `hg_topic_quality()`, `hg_topic_sizes()` | Topic coherence (UMass, NPMI), exclusivity (FREX) and topic sizes | Mimno et al. (2011); Lau, Newman & Baldwin (2014); Roberts, Stewart & Tingley (2019) |
| `hg_stability()`, `hg_agreement()`, `hg_membership()` | Subsampling stability of clusters, agreement between partitions (ARI, AMI, NMI) and soft memberships | Hennig (2007); Hubert & Arabie (1985); Vinh et al. (2010); Bezdek (1981) |
| `hg_cocluster()` | Spectral co-clustering of documents and words | Dhillon (2001) |
| `hg_classify()`, `hg_seeds()` | Transductive classification from a few labels, with class-mass normalization, and seed selection | Zhou, Huang & Schölkopf (2006); Zhu, Ghahramani & Lafferty (2003) |
| `topic_network()` | Network of topics that share documents | Roberts, Stewart & Tingley (2019); Cassi et al. (2017) |
| `hg_sequences()` | Per-actor sequences of clusters from a clustered corpus, for the memory family | Xu, Wickramarathne & Chawla (2016) |

## Naming

Constructors are nouns: `hon()`, `mogen()`, `simplicial()`, `hypergraph()`,
`text_hypergraph()` and the `<source>_hypergraph()` functions. Where a family
has several construction routes, `type =` selects one. All other functions
begin with `hg_`. A measure returns a data frame with one row per node or per
structure, and an inference function returns a result object with estimates,
intervals and p-values. `hg_centrality()`, `hg_communities()` and `hypa()`
take a memory network or a hypergraph.

The tables of every result are read with `hg_get(x, what = )`, where `what`
names a secondary table and the named arguments `sort_by` and `top` order and
shorten it:

```r
library(hypergraphs)

net <- hon(human_long, action = "code", actor = "session_id",
           time = "timestamp", max_order = 3)
hg_get(net)                                # the rules
hg_get(net, what = "nodes")                # the higher-order nodes
hg_get(net, what = "pathways")             # the higher-order pathways

mg <- mogen(human_long, action = "code", actor = "session_id",
            time = "timestamp", max_order = 3)
hg_get(mg, what = "paths", k = 3)          # counts of paths of length 3

hg <- group_hypergraph(human_long, node = "code", hyperedge = "session_id")
hg_get(hg, what = "memberships")           # one row per node in a hyperedge
```

Every function that takes sequences accepts a data frame in long format with
one row per event (`action =`, `actor =`, `time =`, `session =`), a wide data
frame with one sequence per row, a list of sequences, or a fitted network
object that carries its sequences.

## Links between families

The same coded sessions can be read as a memory network, a simplicial complex
and a hypergraph:

```r
net <- hon(human_long, action = "code", actor = "session_id",
           time = "timestamp", max_order = 3)
hg_centrality(net, sort_by = "pagerank")

sc <- simplicial(net, type = "pathway")
hg_betti(sc)

hg <- group_hypergraph(human_long, node = "code", hyperedge = "session_id")
hg_centrality(hg, type = "pagerank")

thg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
                       weight = "tfidf", stop_words = stop_words_en(),
                       min_count = 3L)
topics <- hg_cluster(thg, k = 4, type = "random_walk", seed = 1)
hg_keywords(thg, topics, n = 5, collapse = TRUE)
hg_keywords(thg, topics, n = 5, type = "ctfidf", collapse = TRUE)
```

A text hypergraph is a hypergraph, so every hypergraph function accepts it.
`pairwise_network()` projects a hypergraph to a pairwise network, and
`hg_sequences()` turns a clustered corpus into sequences for the memory
family.

## Data

| Dataset | Content |
|---|---|
| `tutoring_events` | Problem steps of learners with a tutor, one row per event |
| `human_long`, `ai_long` | Coded actions of human-AI pair-programming sessions |
| `covid_abstracts`, `covid_embeddings`, `covid_sample` | COVID-19 education research abstracts (165), their sentence embeddings, and a random sample of 1,000 abstracts |
| `icsid_tribunals` | ICSID arbitration tribunals, one row per seat |
| `gfcc_decisions`, `gfcc_citations` | Decisions of the German Federal Constitutional Court and their citation blocks |
| `forum_posts` | Simulated forum posts whose topics follow second-order memory |
| `ring_sequences`, `ring_communities` | Simulated sequences with planted memory modules, and the planted modules |

## Documentation

The package ships one vignette, `vignette("hypergraphs")`. The other
walkthroughs are articles on the
[package website](https://pak.dynasite.org/hypergraphs/):
[hypergraph analysis of a text corpus](https://pak.dynasite.org/hypergraphs/articles/text-hypergraphs.html), the
[constructions of a text hypergraph](https://pak.dynasite.org/hypergraphs/articles/text-constructions.html),
[topic structure](https://pak.dynasite.org/hypergraphs/articles/covid-topics.html) and
[mixed-membership topics](https://pak.dynasite.org/hypergraphs/articles/topic-mixtures.html) of the COVID-19
education literature, [legal hypergraphs](https://pak.dynasite.org/hypergraphs/articles/legal-hypergraphs.html),
the [Argentina tribunals example](https://pak.dynasite.org/hypergraphs/articles/argentina-tribunals.html) and the
[benchmarks](https://pak.dynasite.org/hypergraphs/articles/benchmarks.html).

## Related packages

hypergraphs covers higher-order structure in a family of network packages
that also includes [psychnets](https://github.com/mohsaqr/psychnets)
(cross-sectional psychometric networks) and
[idiographic](https://github.com/mohsaqr/idiographic) (person-specific
temporal models).

## Installation

hypergraphs is not yet on CRAN. Install it from r-universe, which serves
pre-built binaries:

```r
install.packages("hypergraphs",
                 repos = c("https://mohsaqr.r-universe.dev", "https://cloud.r-project.org"))
```

## References

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M., Patania, A.,
Young, J.-G., & Petri, G. (2020). Networks beyond pairwise interactions:
Structure and dynamics. *Physics Reports*, 874, 1-92.
<https://doi.org/10.1016/j.physrep.2020.05.004>

Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University Press.
<https://doi.org/10.1017/9781108770996>

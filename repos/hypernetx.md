# HyperNetX (HNX) — PNNL (Python)

- **What**: hypergraph analysis library from Pacific Northwest National
  Laboratory; data structures on sparse incidence matrices plus an algorithms
  layer.
- **Verified 2026-08-25** (GitHub source, `hypernetx/algorithms/`): modules =
  `clustering` (`hypergraph_modularity.py` — Kumar-style modularity;
  `laplacians_clustering.py`), `homology`, `generation`, `matching`,
  `metrics`, `temporal`, `concepts`, `embeddings`.
- **Verified 2026-09-29 on 2.4.3** (current on PyPI, released 2026-07-23;
  wheel source inspected): `metrics/s_centrality_measures.py` has
  `s_betweenness_centrality`, `s_closeness_centrality`,
  `s_harmonic_centrality`, `s_harmonic_closeness_centrality`,
  `s_eccentricity`; `clustering/hypergraph_modularity.py` has `modularity`,
  `kumar` (Kumar et al. 2020), `last_step`, `two_section`;
  `generation/generative_models.py` has `erdos_renyi_hypergraph`,
  `chung_lu_hypergraph`, `dcsbm_hypergraph`; `temporal/contagion.py` has
  `discrete_SIR`, `discrete_SIS`, `Gillespie_SIR`; `hif.py` has `to_hif` /
  `from_hif`; `homology` has `betti_numbers`, `homology_basis`. Planned
  oracles for ROADMAP.md Phases 2 and 4 (s-centralities, Kumar modularity,
  generators, HIF).
- **The key fact**: `laplacians_clustering.py` implements the **Hayashi,
  Aksoy, Park & Park (2020) EDVW pipeline exactly** — `prob_trans()`
  (edge-dependent-vertex-weight random-walk transition matrix), `get_pi()`
  (stationary distribution), `norm_lap()` (normalized Laplacian),
  `spec_clus()` (spectral clustering) — citing the paper in the docstrings.
  Sinan Aksoy is at PNNL, so this is the author-adjacent reference
  implementation.
- **Role for us**: primary local oracle for the shipped weighted Laplacian,
  RDC-Spec clustering and EDVW PageRank paths. Its homology module remains a
  possible second TDA oracle. HyperNetX does not supply the Hayashi
  SymNMF/JointNMF branches; hypergraphs now implements RDC-Sym (Algorithm 2),
  J-NMF (Eq. 18), and JS-NMF (Eq. 19), verified directly against their
  paper equations and objective invariants.
- **Links**: https://github.com/pnnl/HyperNetX ·
  https://hypernetx.readthedocs.io · `pip install hypernetx`
- **Full read 2026-10-05 (2.4.3, installed in `~/.virtualenvs/r-reticulate`)**:
  15.4k lines; `classes/hypergraph.py` (2.8k) is the core. Modules not in
  the notes above: `matching` (greedy, iterated-sampling and HEDCS
  approximate matching on d-uniform hypergraphs), `concepts` (formal
  concept lattices, `HypergraphLattice`), `temporal/temporal_paths.py`
  (temporal incidence graph, `temporal_shortest_path()`, temporal line
  graph on edge-ordered hypergraphs), `homology/oat_accelerator.py`
  (rational-coefficient homology through OAT), `reports/descriptive_stats.py`
  (`degree_dist`, `edge_size_dist`, `comp_dist`, `toplex_dist`,
  `s_node_diameter_dist`, `info_dict`), `drawing` (rubber-band, two-column,
  UpSet incidence, storyline). `embeddings/` is an empty package. Class
  methods include `s_connected_components`, `node/edge_diameters`,
  `distance`/`edge_distance` (s-walks), `toplexes`, `collapse_edges/nodes`,
  `equivalence_classes`, set algebra (`sum`, `union`, `intersection`,
  `difference`).
- **Overlap with hypergraphs**: already used as oracle for `hg_laplacian()`,
  `hg_cluster()` spectral, `hg_pagerank()`, `hg_modularity()` and
  `hg_communities(type = "irmm")`. `s_betweenness_centrality()` equals
  `hg_edge_centrality()`; `s_closeness_centrality()` is per-component (see
  LEARNINGS 2026-10-05). Generators cite Aksoy et al. 2017
  (doi:10.1093/comnet/cnx001; Erdos-Renyi, Chung-Lu) and Larremore et al.
  2014 (DCSBM); `random_hypergraph()` has no Chung-Lu or DCSBM type. No
  hypergraphs counterpart for: s-walk distances/diameters, s-components as
  a table, toplexes, edge/node collapsing, contagion (SIR/SIS), temporal
  shortest paths, matching, concept lattices.

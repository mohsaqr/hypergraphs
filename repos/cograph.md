# cograph — graph and plotting engine

- **Relationship**: explicit hypergraphs Import. hypergraphs owns higher-order
  construction and projection; cograph begins when the result is an ordinary
  weighted graph.
- **Used for**: `splot()`, betweenness and closeness on s-line graphs,
  repeated Infomap, community quality and graph rendering.
- **Boundary**: cograph does not construct hypergraphs and hypergraphs does not
  duplicate its ordinary graph algorithms.
- **Collision audit 2026-09-02**: no exported name overlaps with hypergraphs.
- **Version checked locally**: 2.4.5.


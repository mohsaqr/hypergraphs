#!/bin/sh
# Runs every dataset and variant of the HyperGAT reproduction in sequence.
# Usage (from the package root): nohup sh benchmarks/run_hypergat_reproduction_all.sh &
for ds in R8 R52; do
  Rscript benchmarks/run_hypergat_reproduction.R "$ds" lda_official 1:3 || echo "FAILED $ds lda_official"
  Rscript benchmarks/run_hypergat_reproduction.R "$ds" lda 1:10 || echo "FAILED $ds lda"
  Rscript benchmarks/run_hypergat_reproduction.R "$ds" none 1:10 || echo "FAILED $ds none"
done
echo "ALL DONE"

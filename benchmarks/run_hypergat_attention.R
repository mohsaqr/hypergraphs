# Attention weights of HyperGAT on R8 (topics fitted on the training
# documents, seed 1), the counterpart of Figure 5 of Ding et al. (2020).
# Writes the per-(document, word) attention table and the predictions to
# benchmarks/results/hypergat_r8_attention.rds.
# Usage (from the package root): Rscript benchmarks/run_hypergat_attention.R
suppressPackageStartupMessages(library(hypergraphs))
dir <- file.path("local_testing_and_equivalence", "HyperGAT_TextClassification",
                 "data")
docs <- enc2utf8(readLines(file.path(dir, "R8_corpus.txt"), encoding = "latin1",
                           warn = FALSE))
meta <- utils::read.delim(file.path(dir, "R8_labels.txt"), header = FALSE,
                          col.names = c("row", "split", "label"),
                          colClasses = "character", quote = "")
names(docs) <- sprintf("R8_%05d", seq_along(docs))
is_train <- grepl("train", meta$split, fixed = TRUE)
labels <- stats::setNames(meta$label[is_train], names(docs)[is_train])
toks <- table(unlist(hypergraphs:::.thg_tokenize(docs, TRUE)))
stops <- union(stop_words_en(), names(head(sort(toks, decreasing = TRUE), 9)))
fit <- hg_hypergat(docs, labels = labels, stop_words = stops,
                    min_count = 5L, semantic = "lda", seed = 1L)
attention <- hg_get(fit, what = "attention")
predictions <- hg_get(fit)
saveRDS(list(attention = attention, predictions = predictions,
             split = stats::setNames(meta$split, names(docs)),
             label = stats::setNames(meta$label, names(docs)),
             text = docs),
        file.path("benchmarks", "results", "hypergat_r8_attention.rds"))
message("attention written")

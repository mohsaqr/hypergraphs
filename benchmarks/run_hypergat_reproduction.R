# Reproduction of Ding et al. (2020), HyperGAT, Tables 2 and 4 (HyperGAT
# rows) on the official repository's R8 and R52 corpora, which keep the
# sentence structure the sequential hyperedges need. One row per seed is
# written to benchmarks/results/hypergat_reproduction.csv.
#
# Usage: Rscript benchmarks/run_hypergat_reproduction.R <dataset> <variant> <seeds>
#   dataset: R8 | R52
#   variant: lda (native LDA on the training documents, the paper's
#            protocol) | lda_official (the official pickle's topic keywords,
#            fitted on all documents by generate_lda.py) | none (w/o semantic)
#   seeds:   e.g. 1:10
suppressPackageStartupMessages(library(hypergraphs))
source(file.path("benchmarks", "harness.R"))

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[[1]]
variant <- args[[2]]
seeds <- eval(parse(text = args[[3]]))
stopifnot(
  "dataset must be R8 or R52" = dataset %in% c("R8", "R52"),
  "variant must be lda, lda_official or none" =
    variant %in% c("lda", "lda_official", "none")
)

dir <- file.path("local_testing_and_equivalence", "HyperGAT_TextClassification",
                 "data")
# the official preprocess.py decodes the corpus as latin1 (R8 line 5536
# holds the byte 0xfc), so the corpus is read the same way
docs <- enc2utf8(readLines(file.path(dir, paste0(dataset, "_corpus.txt")),
                           encoding = "latin1", warn = FALSE))
meta <- utils::read.delim(file.path(dir, paste0(dataset, "_labels.txt")),
                          header = FALSE, col.names = c("row", "split", "label"),
                          colClasses = "character", quote = "")
stopifnot("corpus and labels must align" = length(docs) == nrow(meta))
names(docs) <- sprintf("%s_%05d", dataset, seq_along(docs))
is_train <- grepl("train", meta$split, fixed = TRUE)
labels <- stats::setNames(meta$label[is_train], names(docs)[is_train])

# the official highbar rule: drop the 9 most frequent tokens
toks <- table(unlist(hypergraphs:::.thg_tokenize(docs, TRUE)))
stops <- union(stop_words_en(), names(head(sort(toks, decreasing = TRUE), 9)))

official_keywords <- if (identical(variant, "lda_official")) {
  kw <- jsonlite::fromJSON(file.path("benchmarks", "results",
                                     paste0(dataset, "_official_lda_keywords.json")))
  unname(kw[order(as.integer(names(kw)))])
}

one_seed <- function(s) {
  semantic_args <- switch(variant,
    lda = list(semantic = "lda"),
    lda_official = list(semantic = "lda", lda_keywords = official_keywords),
    none = list(semantic = "none"))
  t_fit <- system.time(
    fit <- do.call(hg_hypergat, c(list(docs, labels = labels, stop_words = stops,
                                       min_count = 5L, seed = s), semantic_args))
  )[["elapsed"]]
  predicted <- fit$predicted[match(names(docs), fit$node)]
  score <- .bench_score(meta$label[!is_train], predicted[!is_train])
  row <- data.frame(dataset = dataset, variant = variant, seed = s,
                    n_train = sum(is_train), n_test = sum(!is_train),
                    accuracy = score$accuracy, macro_f1 = score$macro_f1,
                    fit_s = t_fit)
  out <- file.path("benchmarks", "results", "hypergat_reproduction.csv")
  utils::write.table(row, out, sep = ",", row.names = FALSE,
                     col.names = !file.exists(out), append = file.exists(out))
  message(sprintf("%s %s seed %d: accuracy %.4f (%.0f s)", dataset, variant, s,
                  score$accuracy, t_fit))
  row
}
invisible(lapply(seeds, one_seed))

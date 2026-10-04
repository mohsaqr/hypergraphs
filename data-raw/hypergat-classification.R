# Regenerate the reference run bundled with the HyperGAT vignette.
# Run from the package root with R8_corpus.txt and R8_labels.txt in data_dir.
# The files are distributed by the official HyperGAT implementation at commit
# 9af29b6837b95344bbc5a567e0673b41c616c88f. Training requires installed torch.
library(hypergraphs)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript data-raw/hypergat-classification.R <R8 data directory>")
}
data_dir <- args[[1L]]
out_dir <- file.path("inst", "extdata", "hypergat-classification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

text <- enc2utf8(readLines(con = file.path(data_dir, "R8_corpus.txt"),
                          encoding = "latin1", warn = FALSE))
meta <- read.delim(file = file.path(data_dir, "R8_labels.txt"), header = FALSE,
                   col.names = c("row", "split", "label"),
                   colClasses = "character", quote = "")
stopifnot(length(text) == nrow(meta))
articles <- data.frame(node = sprintf("R8_%05d", seq_along(text)),
                        text = text, split = meta$split, label = meta$label)
training <- subset(x = articles, subset = split == "train",
                   select = c(node, text, label))
testing <- subset(x = articles, subset = split == "test",
                  select = c(node, text, label))
subjects <- sort(unique(articles$label))
counts <- data.frame(subject = subjects,
  training = as.integer(table(factor(training$label, levels = subjects))),
  test = as.integer(table(factor(testing$label, levels = subjects))))

fit <- hg_hypergat(x = training, column = "text", id = "node",
                   labels = "label", min_count = 5, seed = 1)
evaluation <- predict(object = fit, newdata = testing, labels = "label")
first_article <- "R8_01282"
reference <- list(
  source_commit = "9af29b6837b95344bbc5a567e0673b41c616c88f",
  settings = list(seed = 1L, min_count = 5L, semantic = FALSE),
  counts = counts,
  label_output = capture.output(xtabs(formula = ~ label + split, data = articles)),
  accuracy = hg_get(x = evaluation, what = "accuracy"),
  classes = hg_get(x = evaluation, what = "classes"),
  errors = hg_get(x = evaluation, what = "documents", correct = FALSE, top = 1),
  diagnostics = hg_get(x = fit, what = "hyperedges", node = first_article),
  majority = names(which.max(table(training$label))),
  baseline = mean(testing$label == names(which.max(table(training$label)))),
  scored_n = sum(!is.na(evaluation$predicted)),
  first_article = first_article,
  fit_output = capture.output(print(x = fit)),
  evaluation_output = capture.output(print(x = evaluation))
)
saveRDS(object = reference, file = file.path(out_dir, "results.rds"), version = 2)

save_plot <- function(filename, width, height, draw) {
  grDevices::png(filename = file.path(out_dir, filename),
                width = width, height = height, units = "in", res = 96)
  on.exit(grDevices::dev.off())
  print(draw())
}
save_plot("history.png", 8, 5, function() plot(x = fit, type = "history"))
save_plot("confusion.png", 8, 5, function() plot(x = evaluation))
save_plot("hyperedges.png", 8, 7, function() {
  plot(x = fit, type = "hyperedges", node = first_article, edge_labels = TRUE,
       label_size = 3, titles = FALSE)
})
save_plot("attention.png", 8, 4, function() {
  plot(x = fit, type = "attention", node = first_article)
})

testthat::skip_on_cran()

# ---- hon(group =): the group model of memory networks ----

test_that("each group's network is hon() on that group's sequences", {
  set.seed(1)
  seqs <- replicate(20, sample(c("a", "b", "c", "d"), 12, replace = TRUE),
                    simplify = FALSE)
  labels <- rep(c("x", "y"), 10)
  by_group <- hon(seqs, group = labels, max_order = 2L)
  expect_s3_class(by_group, "hypergraphs_memory_group")
  expect_identical(hg_get(by_group[["x"]]),
                   hg_get(hon(seqs[labels == "x"], max_order = 2L)))
  wide <- as.data.frame(do.call(rbind, seqs), stringsAsFactors = FALSE)
  wide$cohort <- labels
  expect_identical(hg_get(hon(wide, group = "cohort", max_order = 2L)),
                   hg_get(by_group))
  stacked <- hg_get(by_group, what = "nodes")
  expect_identical(names(stacked)[1L], "group")
})

test_that("a long table takes one group label per row, and says so", {
  long <- data.frame(actor = rep(1:4, each = 3),
                     action = rep(c("a", "b", "c"), 4))
  by_row <- hon(long, action = "action", actor = "actor",
                group = rep(c("x", "y"), each = 6), max_order = 1L)
  expect_identical(names(by_row), c("x", "y"))
  err <- tryCatch(hon(long, action = "action", actor = "actor",
                      group = c("x", "x", "y", "y")), error = identity)
  expect_s3_class(err, "hypergraphs_bad_input")
  expect_match(conditionMessage(err), "per row of `data` (12 here)",
               fixed = TRUE)
})

test_that("bad groupings are refused by class", {
  seqs <- list(c("a", "b", "a"), c("b", "a", "b"), c("a", "a", "b"))
  expect_error(hon(seqs, group = c("x", "y")), class = "hypergraphs_bad_input")
  expect_error(hon(seqs, group = rep("x", 3)), class = "hypergraphs_bad_input")
  expect_error(hon(seqs, group = c(NA, "x", "y")),
               class = "hypergraphs_bad_input")
})

test_that("plot() draws one group of a grouped bootstrap", {
  set.seed(1)
  seqs <- replicate(12, sample(c("a", "b", "c"), 10, replace = TRUE),
                    simplify = FALSE)
  by_group <- hon(seqs, group = rep(c("x", "y"), 6), max_order = 2L)
  boot <- hg_bootstrap(by_group, n_boot = 10L, seed = 1L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(boot, group = "y", top = 5L)
  expect_s3_class(p, "ggplot")
  expect_lte(nrow(p$data), 5L)
  expect_error(plot(boot), class = "hypergraphs_bad_input")
  expect_error(plot(boot, group = "zz"), class = "hypergraphs_bad_input")
})

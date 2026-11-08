testthat::skip_on_cran()

.motif_hg <- function(edges) {
  group_hypergraph(
    do.call(rbind, lapply(seq_along(edges), function(i) {
      data.frame(member = edges[[i]], event = paste0("e", i))
    })), "member", "event"
  )
}

# The random comparison with the pairwise loop lives in
# local_testing_and_equivalence/test-equiv-reference-random-2.R.
test_that("motif census counts hand-built Y, T and O sets", {
  # Y: {1,2,3},{1,2,4} share two nodes on {1,2,3,4} -> 1 overlapping pair;
  #    {20,21,22},{20,21,23} likewise -> a second Y.
  # T: three of the four 3-subsets of {5,6,7,8} -> choose(3, 2) = 3 pairs.
  # O: all four 3-subsets of {9,10,11,12} -> choose(4, 2) = 6 pairs.
  # {1,30,31} shares one node with {1,2,3}: no motif. Census Y=2, T=1, O=1.
  edges <- list(c(1L, 2L, 3L), c(1L, 2L, 4L),
                c(20L, 21L, 22L), c(20L, 21L, 23L),
                c(5L, 6L, 7L), c(5L, 6L, 8L), c(5L, 7L, 8L),
                c(9L, 10L, 11L), c(9L, 10L, 12L), c(9L, 11L, 12L),
                c(10L, 11L, 12L), c(1L, 30L, 31L))
  expect_identical(hypergraphs:::.thg_yto_counts(edges),
                   c(Y = 2L, T = 1L, O = 1L))
})

test_that("motif tests carry their draws and plot", {
  h <- .motif_hg(list(c(1, 2, 3), c(1, 2, 4), c(1, 3, 5), c(2, 4, 5)))
  out <- hg_motifs(h, n = 9, seed = 7)
  expect_s3_class(out, "hypergraphs_motifs")
  expect_s3_class(out, "data.frame")
  test <- hg_get(out)
  expect_identical(class(test), "data.frame")
  expect_identical(test$motif, c("Y", "T", "O"))
  draws <- hg_get(out, what = "draws")
  expect_identical(names(draws), c("run", "motif", "count"))
  expect_identical(nrow(draws), 27L)
  expect_equal(mean(draws$count[draws$motif == "Y"]), test$expected[[1L]])
  motif_plot <- plot(out)
  expect_s3_class(motif_plot, "ggplot")
  t_plot <- plot(out, motif = "T")
  expect_s3_class(t_plot, "ggplot")
  counts <- hg_motifs(h, what = "counts")
  expect_identical(class(counts), "data.frame")
})

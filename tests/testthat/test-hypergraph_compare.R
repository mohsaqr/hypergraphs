testthat::skip_on_cran()

# ---- hg_compare() on a hypergraph ------------------------------------------
#
# Node centralities of two groups of hyperedges: the observed values against
# hand values and against hypergraphs built from each group's rows, the
# permutation scheme for each nesting, the bootstrap, and the result class.

# Messages a1-a2 in group "x" hold {p, q}; b1-b2 in group "y" hold {p, r} and
# {q, r}.
.cmp_small <- function() {
  data.frame(
    message = c("a1", "a1", "a2", "a2", "b1", "b1", "b2", "b2"),
    code = c("p", "q", "p", "q", "p", "r", "q", "r"),
    side = rep(c("x", "y"), each = 4),
    stringsAsFactors = FALSE
  )
}

# 24 messages of 6 students, 4 each, in prompts of 2 messages; the level is
# constant within a prompt and varies within a student
.cmp_nested <- function() {
  codes <- list(c("plan", "ask"), c("plan", "check"), c("ask", "check"),
                c("plan", "ask", "check"))
  message <- paste0("m", seq_len(24))
  members <- codes[rep_len(seq_along(codes), 24)]
  data.frame(
    message = rep(message, lengths(members)),
    code = unlist(members),
    student = rep(paste0("s", rep(seq_len(6), each = 4)), lengths(members)),
    prompt = rep(paste0("p", rep(seq_len(12), each = 2)), lengths(members)),
    level = rep(rep(c("high", "high", "low", "low"), 6), lengths(members)),
    stringsAsFactors = FALSE
  )
}

test_that("degree differences equal the hand values", {
  hg <- hypergraph(.cmp_small(), node = "code", hyperedge = "message")
  cmp <- hg_compare(hg, group = "side", type = "degree", n_perm = 19,
                    n_boot = 19, seed = 1)
  tab <- hg_get(cmp)
  expect_identical(tab$node, c("p", "q", "r"))
  # x: p and q in both messages, r in none; y: p, q in one of two, r in both
  expect_equal(tab$x, c(1, 1, 0))
  expect_equal(tab$y, c(0.5, 0.5, 1))
  expect_equal(tab$diff, c(0.5, 0.5, -1))
  expect_identical(tab$higher, c("x", "x", "y"))
  expect_named(tab, c("node", "measure", "x", "y", "diff", "lower", "upper",
                      "p_value", "p_adj", "significant", "higher"))
})

test_that("observed centralities equal those of each group's own hypergraph", {
  d <- .cmp_nested()
  hg <- hypergraph(d, node = "code", hyperedge = "message")
  cmp <- hg_compare(hg, group = "level", cluster = "student",
                    n_perm = 9, n_boot = 9, seed = 1)
  tab <- hg_get(cmp)
  lapply(c("high", "low"), \(g) {
    own <- hypergraph(subset(d, level == g), node = "code",
                      hyperedge = "message")
    cen <- hg_centrality(own, type = c("clique", "pagerank"))
    deg <- hg_measures(own)
    lapply(c("clique", "pagerank"), \(m) {
      rows <- subset(tab, measure == m)
      expect_equal(rows[[g]], cen[[m]][match(rows$node, cen$node)])
    })
    rows <- subset(tab, measure == "degree")
    expect_equal(rows[[g]],
                 deg$hyperdegree[match(rows$node, deg$node)] / 12)
  })
})

test_that("swapping the groups negates the differences and keeps the p-values", {
  hg <- hypergraph(.cmp_nested(), node = "code", hyperedge = "message")
  one <- hg_get(hg_compare(hg, group = "level", groups = c("high", "low"),
                           cluster = "student", n_perm = 49, n_boot = 49,
                           seed = 3))
  two <- hg_get(hg_compare(hg, group = "level", groups = c("low", "high"),
                           cluster = "student", n_perm = 49, n_boot = 49,
                           seed = 3))
  expect_equal(two$diff, -one$diff)
  expect_equal(two$lower, -one$upper)
  expect_equal(two$p_value, one$p_value)
  expect_identical(two$higher[one$higher != "equal"],
                   one$higher[one$higher != "equal"])
})

test_that("groups with the same hyperedges differ by zero with p = 1", {
  d <- .cmp_small()
  twin <- rbind(transform(d, side = "x"),
                transform(d, side = "y", message = paste0(d$message, "'")))
  hg <- hypergraph(twin, node = "code", hyperedge = "message")
  tab <- hg_get(hg_compare(hg, group = "side", n_perm = 19, n_boot = 19,
                           seed = 2))
  expect_equal(tab$diff, rep(0, nrow(tab)))
  expect_equal(tab$p_value, rep(1, nrow(tab)))
  expect_identical(unique(tab$higher), "equal")
})

test_that("a planted difference is detected and the seed reproduces it", {
  set.seed(11)
  n <- 40
  members <- c(replicate(n, c("a", "b"), simplify = FALSE),
               replicate(n, sample(c("a", "b", "c", "d"), 2), simplify = FALSE))
  d <- data.frame(message = rep(paste0("m", seq_len(2 * n)), lengths(members)),
                  code = unlist(members),
                  level = rep(rep(c("one", "two"), each = n), lengths(members)))
  hg <- hypergraph(d, node = "code", hyperedge = "message")
  before <- .Random.seed
  cmp <- hg_compare(hg, group = "level", type = "degree", n_perm = 199,
                    n_boot = 199, seed = 5)
  expect_identical(.Random.seed, before)
  tab <- hg_get(cmp)
  expect_true(all(subset(tab, node %in% c("c", "d"))$significant))
  expect_identical(subset(tab, node == "c")$higher, "two")
  again <- hg_compare(hg, group = "level", type = "degree", n_perm = 199,
                      n_boot = 199, seed = 5)
  expect_identical(hg_get(again), tab)
})

test_that("labels are permuted at the level the group varies", {
  d <- .cmp_nested()
  hg <- hypergraph(d, node = "code", hyperedge = "message")
  within <- hg_compare(hg, group = "level", cluster = "student",
                       n_perm = 9, n_boot = 9, seed = 1)
  expect_identical(within$permuted, "units within clusters")
  expect_identical(within$n_clusters, 6L)
  # the level constant within each student: students carry the labels
  between_data <- transform(d, level = ifelse(student %in% c("s1", "s2", "s3"),
                                              "high", "low"))
  between <- hg_compare(hypergraph(between_data, node = "code",
                                   hyperedge = "message"),
                        group = "level", cluster = "student", n_perm = 9,
                        n_boot = 9, seed = 1)
  expect_identical(between$permuted, "clusters")
  # a restricted permutation keeps every cluster's labels: the group sizes of
  # every permuted labelling equal the observed ones
  design <- .hg_cmp_design(hg, "level", NULL, "student", "prompt")
  set.seed(4)
  lapply(seq_len(20), \(i) {
    shuffled <- .hg_cmp_shuffle(design)
    expect_identical(table(shuffled, design$cluster_of_edge, dnn = NULL),
                     table(design$labels, design$cluster_of_edge, dnn = NULL))
    # a prompt keeps one label
    expect_true(all(tapply(shuffled, design$unit_of_edge,
                           \(v) length(unique(v))) == 1L))
  })
})

test_that("a bootstrap draw with one group is left out with a warning", {
  d <- transform(.cmp_small(), student = ifelse(side == "x", "s1", "s2"))
  hg <- hypergraph(d, node = "code", hyperedge = "message")
  expect_warning(
    cmp <- hg_compare(hg, group = "side", cluster = "student", n_perm = 9,
                      n_boot = 40, seed = 1),
    class = "hypergraphs_undefined_statistic")
  expect_lt(cmp$n_boot_used, 40L)
  expect_identical(cmp$permuted, "clusters")
})

test_that("hg_compare() on a hypergraph refuses bad input with classed errors", {
  d <- .cmp_nested()
  hg <- hypergraph(d, node = "code", hyperedge = "message")
  bad <- list(
    list(group = "nope"), list(group = "level", cluster = "nope"),
    list(group = "level", groups = c("high", "mid")),
    list(group = "level", groups = c("high", "high")),
    list(group = "level", type = "katz"),
    list(group = "level", type = c("degree", "degree")),
    list(group = "level", n_perm = 1), list(group = "level", n_boot = 1.5),
    list(group = "level", level = 1), list(group = "level", alpha = 0),
    list(group = "level", seed = 1.5),
    # the memory-network method's argument
    list(group = "level", parallel = TRUE),
    # the level varies within a student, so a student is not a unit
    list(group = "level", unit = "student"),
    # a prompt holds messages of one student; a student spans prompts
    list(group = "level", cluster = "prompt", unit = "student")
  )
  lapply(bad, \(args) {
    call_args <- utils::modifyList(list(n_perm = 9, n_boot = 9), args)
    expect_error(do.call(hg_compare, c(list(hg), call_args)),
                 class = "hypergraphs_bad_input")
  })
  expect_error(hg_compare(hg, n_perm = 9, n_boot = 9),
               class = "hypergraphs_bad_input")
  expect_error(hg_compare(hg, group = "level", n_perm = 1, n_boot = 9),
               class = "hypergraphs_bad_input")
  three <- hypergraph(transform(d, level = rep_len(c("a", "b", "c"), nrow(d))),
                      node = "code", hyperedge = "message")
  expect_error(hg_compare(three, group = "level", n_perm = 9, n_boot = 9),
               class = "hypergraphs_bad_input")
  expect_error(hg_compare(data.frame(a = 1)), class = "hypergraphs_bad_input")
})

test_that("hyperedges without a group are left out with a warning", {
  d <- transform(.cmp_nested(), level = ifelse(message == "m1", NA, level))
  hg <- hypergraph(d, node = "code", hyperedge = "message")
  expect_warning(
    cmp <- hg_compare(hg, group = "level", n_perm = 9, n_boot = 9, seed = 1),
    class = "hypergraphs_dropped_rows")
  expect_identical(sum(cmp$n_hyperedges), 23L)
})

test_that("the comparison is read, printed, summarised and plotted", {
  hg <- hypergraph(.cmp_nested(), node = "code", hyperedge = "message")
  cmp <- hg_compare(hg, group = "level", cluster = "student",
                    n_perm = 19, n_boot = 19, seed = 1)
  tab <- hg_get(cmp)
  expect_identical(class(tab), "data.frame")
  expect_identical(nrow(tab), 9L)
  expect_identical(unique(hg_get(cmp, measure = "clique")$measure), "clique")
  sorted <- hg_get(cmp, sort_by = "abs_diff")
  expect_identical(abs(sorted$diff), sort(abs(tab$diff), decreasing = TRUE))
  expect_identical(nrow(hg_get(cmp, top = 2)), 2L)
  expect_true(all(hg_get(cmp, significant = TRUE)$significant))
  expect_error(hg_get(cmp, measure = "katz"), class = "hypergraphs_bad_input")
  expect_error(hg_get(cmp, what = "rules"), class = "hypergraphs_bad_input")
  expect_error(hg_get(cmp, sort_by = "node"), class = "hypergraphs_bad_input")
  expect_output(print(cmp), "Hypergraph comparison of `level`: high")
  expect_output(print(cmp), "permutations of units within clusters")
  s <- summary(cmp)
  expect_s3_class(s, "hypergraphs_summary")
  expect_identical(s$by_measure$measure, c("degree", "clique", "pagerank"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(cmp), "ggplot")
})

test_that("hg_compare() dispatches on group models and hypergraphs", {
  seqs <- list(c("a", "b", "c", "a", "b"), c("a", "b", "c", "a", "b"),
               c("x", "b", "d", "x", "b"), c("x", "b", "d", "x", "b"))
  model <- hon(seqs, group = c("g1", "g1", "g2", "g2"), max_order = 2)
  expect_s3_class(hg_compare(model, n_perm = 9, seed = 1),
                  "hypergraphs_comparison")
  expect_error(hg_compare(model, group = "level", n_perm = 9),
               class = "hypergraphs_bad_input")
  hg <- hypergraph(.cmp_nested(), node = "code", hyperedge = "message")
  expect_s3_class(hg_compare(hg, group = "level", n_perm = 9, n_boot = 9,
                             seed = 1), "hypergraphs_hg_comparison")
})

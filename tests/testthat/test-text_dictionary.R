# hg_dictionary() and text_hypergraph(separator =)

papers <- data.frame(
  keywords = c("Medical Education; Online Learning",
               "teacher education;online learning",
               "medical students; assessment; Medical education",
               "online learning; assessment",
               "teacher education; medical students",
               "  Higher   Education ;; assessment "),
  stringsAsFactors = FALSE
)

test_that("separator reads a delimited field as whole-phrase terms", {
  kw <- text_hypergraph(papers, column = "keywords", separator = ";")
  expect_setequal(colnames(kw$incidence),
                  c("medical education", "online learning",
                    "teacher education", "medical students", "assessment",
                    "higher education"))
  # case folded, whitespace collapsed, empty pieces dropped
  expect_identical(sum(kw$incidence[, "medical education"] != 0), 2L)
  expect_identical(sum(kw$incidence[, "higher education"] != 0), 1L)
  # stop_words drops whole terms, min_count rare terms
  without <- text_hypergraph(papers, column = "keywords", separator = ";",
                             stop_words = "online learning")
  expect_false("online learning" %in% colnames(without$incidence))
  common <- suppressMessages(text_hypergraph(papers, column = "keywords",
                                             separator = ";", min_count = 2L))
  expect_false("higher education" %in% colnames(common$incidence))
  # edge punctuation, quotes and no-break spaces go; internal hyphens stay
  messy <- data.frame(k = c("\u2013 COVID-19; \u2018becoming'; K-12.",
                            "covid-19;\u00a0k-12 ; e-learning"))
  terms <- text_hypergraph(messy, column = "k", separator = ";")
  expect_setequal(colnames(terms$incidence),
                  c("covid-19", "becoming", "k-12", "e-learning"))
  expect_identical(hg_get(terms, what = "vocabulary", sort_by = "count",
                          top = 2)$word, c("covid-19", "k-12"))
  expect_error(text_hypergraph(papers, column = "keywords", separator = ";",
                               construction = "window"),
               class = "hypergraphs_bad_input")
  expect_error(text_hypergraph(papers, column = "keywords", separator = ""),
               class = "hypergraphs_bad_input")
})

test_that("hg_dictionary() labels each node by its most matched category", {
  kw <- text_hypergraph(papers, column = "keywords", separator = ";")
  dictionary <- list(health = c("medical education", "medical students"),
                     teachers = "teacher education")
  expect_message(coded <- hg_dictionary(kw, dictionary),
                 "3 of 6 nodes labelled; 2 match no term and 1 tie")
  expect_identical(coded$node, c("doc_1", "doc_3", "doc_2"))
  expect_identical(coded$label, c("health", "health", "teachers"))
  # doc_3 carries both health terms; doc_5 ties one health and one teacher
  expect_identical(coded$n_terms, c(1L, 2L, 1L))
  expect_identical(coded$terms[[2]], "medical education; medical students")
  expect_false("doc_5" %in% coded$node)
  # the table is a labels input as it stands
  expect_identical(.thg_labels_input(coded),
                   c(doc_1 = "health", doc_3 = "health", doc_2 = "teachers"))
})

test_that("hg_dictionary() refuses bad dictionaries and warns on unknown terms", {
  kw <- text_hypergraph(papers, column = "keywords", separator = ";")
  expect_error(hg_dictionary(kw, list("medical education")),
               class = "hypergraphs_bad_input")
  expect_error(hg_dictionary(kw, list(a = "assessment", b = "assessment")),
               class = "hypergraphs_bad_input")
  expect_error(suppressWarnings(hg_dictionary(kw, list(a = "nothing"))),
               class = "hypergraphs_bad_input")
  expect_warning(suppressMessages(hg_dictionary(
    kw, list(health = c("medical education", "nursing education")))),
    class = "hypergraphs_unmatched_terms")
})

test_that("unlabelled documents are their own split", {
  corpus <- data.frame(
    text = c("soup salt onion broth", "salt soup onion", "broth soup salt",
             "onion salt broth", "stars sky moon night", "sky stars moon",
             "night sky stars", "moon night sky", "soup night", "stars soup"),
    theme = c(rep(c("food", "sky"), each = 4), NA, NA)
  )
  hg <- text_hypergraph(corpus, column = "text")
  fit <- hg_classify(hg, labels = "theme", holdout = 0.5)
  unlabelled <- hg_get(fit, split = "unlabelled")
  expect_identical(unlabelled$node, c("doc_10", "doc_9"))
  expect_true(all(is.na(unlabelled$label)))
  expect_identical(nrow(hg_get(fit, split = "test")) +
                     nrow(hg_get(fit, split = "train")) + 2L, 10L)
})

test_that("hg_get() reads the documents named by a table, in its order", {
  corpus <- data.frame(text = c("soup salt", "stars moon", "soup moon"),
                       title = c("Food", "Sky", "Both"))
  hg <- text_hypergraph(corpus, column = "text")
  picked <- data.frame(node = c("doc_3", "doc_1"), predicted = c("x", "y"))
  out <- hg_get(hg, what = "documents", node = picked)
  expect_identical(out$doc, c("doc_3", "doc_1"))
  expect_identical(out$title, c("Both", "Food"))
  expect_error(hg_get(hg, node = "doc_1"), class = "hypergraphs_bad_input")
})

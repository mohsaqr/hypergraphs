skip_on_cran()

# non-ASCII input is written as \u escapes so this file is portable:
# \u2019 curly apostrophe, \u00a9 copyright sign, \u00e2\u20ac\u2122 the
# three-character mojibake of a curly apostrophe
messy <- c(
  list = "1. The programme raised attainment in 2020, see section 3.2.1.",
  journal = "(3) OJ No L 297, 24.11.1979, p. 1.",
  html = "Peer&nbsp;support <b>improved</b> outcomes [12] across schools.",
  url = "Full report at https://example.org/study.pdf and doi:10.1/x.",
  notice = paste0("Students\u2019 views shifted. \u00a9 2021 Informa UK ",
                  "Limited, trading as Taylor & Francis Group."),
  mojibake = "Teachers\u00e2\u20ac\u2122 workload rose by 25% [3, 4].",
  entity = "Learning &amp; teaching &#8217;online&#8217; in 2021.",
  rights = "Simulation training works. All rights reserved."
)

test_that("clean_text repairs and strips each class of non-content", {
  out <- clean_text(messy)
  expect_identical(length(out), length(messy))
  expect_identical(names(out), names(messy))
  expect_identical(out[["list"]],
                   "The programme raised attainment in, see section.")
  expect_identical(out[["html"]],
                   "Peer support improved outcomes across schools.")
  expect_identical(out[["url"]], "Full report at and.")
  expect_identical(out[["notice"]], "Students' views shifted.")
  expect_identical(out[["mojibake"]], "Teachers' workload rose by.")
  expect_identical(out[["entity"]], "Learning & teaching 'online' in.")
  expect_identical(out[["rights"]], "Simulation training works.")
  # every result is plain ASCII
  expect_true(all(vapply(out, \(s) all(utf8ToInt(s) < 128), logical(1))))
})

test_that("clean_text switches, stop words, min_content and NA behave", {
  keep_numbers <- clean_text(messy[["list"]], numbers = FALSE)
  expect_true(grepl("2020", keep_numbers, fixed = TRUE))
  expect_false(grepl("3.2.1", keep_numbers, fixed = TRUE))
  raw_html <- clean_text(messy[["html"]], html = FALSE)
  expect_true(grepl("<b>", raw_html, fixed = TRUE))
  kept_copyright <- clean_text("Copyright 2020 Elsevier. Great work.",
                               copyright = FALSE)
  expect_identical(kept_copyright, "Copyright Elsevier. Great work.")
  leading_notice <- clean_text("Copyright 2020 Elsevier. Great work.")
  expect_identical(leading_notice, "")
  trailing_notice <- clean_text("Great work. Copyright 2020 Elsevier.")
  expect_identical(trailing_notice, "Great work.")
  sw <- clean_text("The cat sat on the mat", stop_words = c("the", "on"))
  expect_identical(sw, "cat sat mat")
  # the content floor: the journal citation is nearly all digits/punctuation
  floored <- clean_text(messy, min_content = 0.5)
  expect_identical(unname(floored[["journal"]]), "")
  expect_true(nzchar(floored[["html"]]))
  with_na <- clean_text(c("a", NA_character_))
  expect_identical(with_na, c("a", ""))
  empty <- clean_text(character(0))
  expect_identical(empty, character(0))
  # idempotent: cleaning clean text changes nothing
  once <- clean_text(messy)
  twice <- clean_text(once)
  expect_identical(twice, once)
})

test_that("clean_text keeps data.frame rows aligned and feeds text_hypergraph", {
  df <- data.frame(id = c("a", "b", "c"),
                   text = c(messy[["html"]], messy[["journal"]],
                            messy[["notice"]]),
                   year = c(2020L, 1979L, 2021L))
  out <- clean_text(df, column = "text", min_content = 0.5)
  expect_identical(dim(out), dim(df))
  expect_identical(out$id, df$id)
  expect_identical(out$year, df$year)
  expect_identical(out$text[[2L]], "")
  # the empty row is dropped by the constructor, with its named warning
  expect_warning(hg <- text_hypergraph(out, column = "text", id = "id"),
                 class = "hypergraphs_dropped_documents")
  expect_setequal(hg$nodes, c("a", "c"))
  expect_error(clean_text(df), class = "hypergraphs_bad_input")
  expect_error(clean_text(df, column = "nope"), class = "hypergraphs_bad_input")
  expect_error(clean_text("x", column = "text"), class = "hypergraphs_bad_input")
  expect_error(clean_text(1), class = "hypergraphs_bad_input")
  expect_error(clean_text("x", html = NA), class = "hypergraphs_bad_input")
  expect_error(clean_text("x", min_content = 2), "min_content")
})

test_that("clean_text `remove` patterns empty placeholders before the floor", {
  x <- c("[No abstract available]", "A real abstract about schools.")
  out <- clean_text(x, remove = "no abstract available")
  expect_identical(out, c("", "A real abstract about schools."))
  # without the pattern the placeholder survives as ordinary text
  kept <- clean_text(x)
  expect_identical(kept[[1L]], "[No abstract available]")
  # several patterns, case-insensitive, regex
  out2 <- clean_text("Published by ELSEVIER Ltd; Data: Table 2.",
                     remove = c("published by \\w+ ltd", "table \\d"))
  expect_identical(out2, "Data:.")
  expect_error(clean_text("x", remove = NA_character_), "remove")
})

# --- min_chars -------------------------------------------------------------

test_that("clean_text(min_chars) drops short words and keeps longer ones", {
  x <- "M. Wathelet and J. N. Cunha Rodrigues wrote on energy"
  expect_identical(clean_text(x, min_chars = 3),
                   "Wathelet and. Cunha Rodrigues wrote energy")
  # default is a no-op: cleaning without min_chars is unchanged
  expect_identical(clean_text(x), clean_text(x, min_chars = 0))
  expect_identical(clean_text(x), clean_text(x, min_chars = 1))
})

test_that("clean_text(min_chars) never splits a word at an apostrophe", {
  expect_identical(clean_text("children's views", min_chars = 3),
                   "children's views")
})

test_that("clean_text(min_chars) is monotone in min_chars", {
  x <- "a bo cat food plates"
  kept <- vapply(1:6, \(k) length(strsplit(clean_text(x, min_chars = k),
                                           " ")[[1]]), integer(1))
  expect_false(is.unsorted(rev(kept)))
})

test_that("clean_text rejects a bad min_chars", {
  expect_error(clean_text("text", min_chars = -1))
  expect_error(clean_text("text", min_chars = c(2, 3)))
})

test_that("clean_text() removes the word Copyright before the sign", {
  sign <- "\u00a9"
  expect_identical(
    clean_text(paste0("Teachers adapted. Copyright ", sign,
                      " 2020 Elsevier Ltd.")),
    "Teachers adapted."
  )
  expect_identical(clean_text("Teachers adapted. Copyright (c) 2021 Wiley."),
                   "Teachers adapted.")
  # a sign without the word is removed as before
  expect_identical(clean_text(paste0("Teachers adapted. ", sign, " 2022 IEEE.")),
                   "Teachers adapted.")
})

test_that("clean_text(boilerplate = TRUE) removes publisher boilerplate anywhere", {
  text <- c(
    "Published by Elsevier Ltd. Teachers moved online during the lockdown.",
    "Students valued feedback. Taylor & Francis Group, LLC, all rights reserved.",
    "Taylor (2019) studied tutors under exclusive licence to Springer Nature.",
    "The American Chemical Society and Division of Chemical Education, Inc."
  )
  out <- clean_text(text, boilerplate = TRUE)
  expect_identical(out[[1]], "Teachers moved online during the lockdown.")
  expect_false(any(grepl("(?i)elsevier|francis|llc|rights|springer|chemical",
                         out, perl = TRUE)))
  # a surname that also names a publisher survives outside the full name
  expect_true(grepl("Taylor", out[[3]], fixed = TRUE))
  expect_true(grepl("studied tutors", out[[3]], fixed = TRUE))
  # off by default: the result equals cleaning without the switch
  expect_identical(clean_text(text), clean_text(text, boilerplate = FALSE))
  expect_true(grepl("Elsevier", clean_text(text)[[1]], fixed = TRUE))
  expect_error(clean_text(text, boilerplate = "yes"),
               class = "hypergraphs_bad_input")
})

test_that("clean_text() repairs UTF-8 read as Windows-1252, emoji included", {
  thumbs <- "\U0001F44D"
  # build the garble from the bytes, as an export does
  bytes <- charToRaw(enc2utf8(thumbs))
  garble <- paste0(vapply(bytes, \(b) {
    out <- iconv(rawToChar(b), from = "CP1252", to = "UTF-8")
    if (is.na(out)) intToUtf8(as.integer(b)) else out
  }, character(1L)), collapse = "")
  expect_false(identical(garble, thumbs))
  expect_identical(clean_text(paste("Nice work", garble, "(2)"),
                              numbers = FALSE),
                   paste("Nice work", thumbs, "(2)"))
  # an accented word that is not mojibake is left as it is
  expect_identical(clean_text("Gu\u00f0r\u00fan wrote the draft"),
                   "Gu\u00f0r\u00fan wrote the draft")
  # the garbled e-acute and curly apostrophe still repair
  expect_identical(clean_text("caf\u00c3\u00a9 it\u00e2\u20ac\u2122s open"),
                   "caf\u00e9 it's open")
})

# ---- audit regressions (2026-10-06) -----------------------------------------

test_that("number removal keeps alphanumeric words whole (TXT-03)", {
  kept <- c("covid19", "p53", "abc123def", "123abc", "covid-19", "sars-cov-2",
            "3.5abc", "5three")
  expect_identical(clean_text(kept, citations = FALSE), kept)
  removed <- clean_text(c("in 2020 we saw 45% of 1,000 students",
                          "the 3rd and 21st cases", "values -5 and +3.5",
                          "pages 10-20 here"), citations = FALSE)
  expect_identical(removed, c("in we saw of students", "the and cases",
                              "values and", "pages here"))
  mixed <- clean_text("covid19 rose 12% in 2021 (p53)", citations = FALSE)
  expect_identical(mixed, "covid19 rose in (p53)")
  # cleaning is idempotent
  again <- clean_text(c(kept, removed, mixed), citations = FALSE)
  expect_identical(clean_text(again, citations = FALSE), again)
})

test_that("numeric entities that name no character decode to U+FFFD (TXT-09)", {
  replacement <- intToUtf8(65533L)
  out <- clean_text(c("safe &#55296; text", "hex &#xDFFF; tail",
                      "valid &#233; and &#x1F600; and &#65;",
                      "mix &#xD800; and &#8217;ok&#8217;"), numbers = FALSE)
  expect_identical(out[[1L]], paste("safe", replacement, "text"))
  expect_identical(out[[2L]], paste("hex", replacement, "tail"))
  expect_identical(out[[3L]], paste("valid", intToUtf8(233L), "and",
                                    intToUtf8(128512L), "and A"))
  # the valid curly quotes decode, then normalise to straight quotes
  expect_identical(out[[4L]], paste0("mix ", replacement, " and 'ok'"))
  expect_false(anyNA(out))
  expect_length(clean_text("safe &#55296; text"), 1L)
})

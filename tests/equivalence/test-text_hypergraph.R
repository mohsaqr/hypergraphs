# text_hypergraph(), stop_words_en(): a corpus as a document-word hypergraph. With weight = "tfidf" a cell is the
# word's count in the document times its smoothed inverse document frequency ln((1 + N) / (1 + df)) + 1, without
# normalisation (scikit-learn's TfidfVectorizer(smooth_idf = True, norm = None)). References: scikit-learn, the
# counts and the formula in base R, and the Snowball English stop-word list as tm and stopwords ship it.

docs <- c(d1 = "the cat sat on the mat", d2 = "the dog sat on the log", d3 = "cats and dogs are pets",
          d4 = "the mat and the log", d5 = "a dog and a cat met on a mat")

test_that("counts and tf-idf weights equal the token counts and the smoothed idf formula", {
  tokens <- strsplit(tolower(docs), "[^a-z]+")
  counted <- do.call(rbind, Map(\(words, doc) {
    table_words <- table(words)
    data.frame(doc = doc, word = names(table_words), count = as.vector(table_words))
  }, tokens, names(docs)))
  frequency <- table(counted$word)
  idf <- log((1 + length(docs)) / (1 + as.vector(frequency[counted$word]))) + 1
  # One logarithm and one product per cell.
  expect_reference(hg_get(text_hypergraph(docs, weight = "tfidf"), what = "weights"),
                   transform(counted, weight = count * idf),
                   reference = "independent: token counts and ln((1 + N) / (1 + df)) + 1, base R",
                   covers = c("text_hypergraph()", "hg_get.text_hypergraph()"), quantity = "tf-idf weights", tolerance = 1e-12,
                   key = c("doc", "word"))
})

test_that("tf-idf weights equal scikit-learn's TfidfVectorizer without normalisation", {
  skip_without_python(c("scikit-learn" = "sklearn"))
  script <- c(
    "import sys, json",
    "from sklearn.feature_extraction.text import TfidfVectorizer",
    "spec = json.load(open(sys.argv[1]))",
    "vectorizer = TfidfVectorizer(smooth_idf=True, norm=None, token_pattern=r'(?u)\\b\\w+\\b')",
    "X = vectorizer.fit_transform(spec['text']).tocoo()",
    "words = vectorizer.get_feature_names_out()",
    "json.dump({'doc': [spec['doc'][i] for i in X.row], 'word': [str(words[j]) for j in X.col], 'weight': [float(v) for v in X.data]},",
    "          open(sys.argv[2], 'w'))")
  theirs <- as.data.frame(ref_python_run(script, list(doc = names(docs), text = unname(docs))))
  expect_reference(hg_get(text_hypergraph(docs, weight = "tfidf"), what = "weights"), theirs, reference = "scikit-learn",
                   covers = c("text_hypergraph()", "hg_get.text_hypergraph()"), quantity = "tf-idf weights", tolerance = 1e-12,
                   key = c("doc", "word"),
                   note = paste(ref_python_version("scikit-learn"), "TfidfVectorizer(smooth_idf = True, norm = None), single-character tokens kept"))
})

test_that("the Snowball stop-word list equals tm's and differs from stopwords' by one word", {
  skip_without_reference(c("tm", "stopwords"))
  ours <- data.frame(word = stop_words_en("snowball"), listed = TRUE)
  expect_reference(ours, data.frame(word = tm::stopwords("english"), listed = TRUE), reference = "tm", covers = "stop_words_en()",
                   quantity = "Snowball English stop words", tolerance = 0, key = "word")
  expect_known_difference(ours, data.frame(word = stopwords::stopwords("en", source = "snowball"), listed = TRUE),
                          reference = "stopwords", covers = "stop_words_en()", quantity = "Snowball English stop words",
                          tolerance = 0, key = "word",
                          reason = "The stopwords package's Snowball list has one word more ('will'); hypergraphs' list is the 174 words of tm::stopwords('english').")
})

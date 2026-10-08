# English function-word stop list

A fixed list of English function words for the `stop_words` argument of
[`text_hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/text_hypergraph.md)
and
[`clean_text()`](https://pak.dynasite.org/hypergraphs/reference/clean_text.md).
`type = "minimal"` (the default) is a small list of articles,
prepositions, conjunctions and auxiliaries, deliberately minimal and
versioned with the package. `type = "snowball"` is the English stop list
of the Snowball project (Porter 2001), 174 words that add the personal
pronouns, the forms of "do" and the contractions ("i", "you", "don't",
"i'm"), which chat messages and other informal text are full of.
Corpus-specific boilerplate (e.g. "study", "results" in an abstract
corpus) is added by the caller, as in
`c(stop_words_en(), "study", "results")`.

## Usage

``` r
stop_words_en(type = c("minimal", "snowball"))
```

## Arguments

- type:

  `"minimal"` (default) or `"snowball"`.

## Value

A sorted character vector of lowercase English function words. Raises
`hypergraphs_bad_input` for an unknown `type`.

## References

Manning, C. D., Raghavan, P., & Schütze, H. (2008). *Introduction to
Information Retrieval*. Cambridge University Press.
[doi:10.1017/CBO9780511809071](https://doi.org/10.1017/CBO9780511809071)

Porter, M. F. (2001). Snowball: A language for stemming algorithms.
<https://snowballstem.org/texts/introduction.html>

## Examples

``` r
hg <- text_hypergraph(
  c(a = "the salt and the soup", b = "the soup and the stars"),
  stop_words = stop_words_en()
)
hg_get(hg, what = "vocabulary")
#>    word count doc_freq
#> 1  salt     1        1
#> 2  soup     2        2
#> 3 stars     1        1
length(stop_words_en(type = "snowball"))
#> [1] 174
```

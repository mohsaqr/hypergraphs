# Label nodes with a dictionary of terms

Codes each node by the hyperedges it belongs to. In a text hypergraph
with documents as nodes, the hyperedges of a document are its words, or
with `text_hypergraph(separator = ";")` its keywords, so a dictionary
that lists the terms of each category labels every document that carries
one of them. This is dictionary-based content analysis (Grimmer &
Stewart 2013): the categories and their terms are the analyst's, and the
labelling is a count. A node takes the category with the most matched
terms. A node that matches no term, or matches two categories equally,
takes no label, and a message reports how many. The result is the
`labels` input of
[`hg_classify()`](https://pak.dynasite.org/hypergraphs/reference/hg_classify.md)
and
[`hg_hypergat()`](https://pak.dynasite.org/hypergraphs/reference/hg_hypergat.md),
which spread or learn the labels to the nodes the dictionary cannot
code.

## Usage

``` r
hg_dictionary(hg, dictionary)
```

## Arguments

- hg:

  A `net_hg`, typically from
  [`text_hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/text_hypergraph.md).

- dictionary:

  A named list of character vectors, one per category, each holding the
  terms of that category as they appear among the hyperedge names
  (lowercase for a text hypergraph built with `lowercase = TRUE`). A
  term may belong to one category only.

## Value

A base `data.frame`, one row per labelled node: `node`, `label` (the
category), `n_terms` (how many of the category's terms the node carries)
and `terms` (those terms, joined by `"; "`), ordered by `label` and
`node`. Raises `hypergraphs_bad_input` for a malformed dictionary and
warns with `hypergraphs_unmatched_terms` when a term of the dictionary
is not a hyperedge of `hg`.

## References

Grimmer, J., & Stewart, B. M. (2013). Text as data: The promise and
pitfalls of automatic content analysis methods for political texts.
*Political Analysis*, 21(3), 267-297.
[doi:10.1093/pan/mps028](https://doi.org/10.1093/pan/mps028)

## Examples

``` r
papers <- data.frame(
  keywords = c("medical education; online learning",
               "teacher education; online learning",
               "medical students; assessment",
               "online learning; assessment")
)
kw <- text_hypergraph(papers, column = "keywords", separator = ";")
hg_dictionary(kw, dictionary = list(
  health = c("medical education", "medical students"),
  teachers = "teacher education"
))
#> 3 of 4 nodes labelled; 1 match no term and 0 tie between categories
#>    node    label n_terms             terms
#> 1 doc_1   health       1 medical education
#> 2 doc_3   health       1  medical students
#> 3 doc_2 teachers       1 teacher education
```

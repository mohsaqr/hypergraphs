# Simulated forum posts whose topics follow second-order memory

A corpus whose sequence of topics can only be described with memory.
Thirty students each write twelve posts. Each post is about one of three
topics, cooking, astronomy or gardening, and holds seven words drawn
with replacement from the vocabulary of its topic and one everyday word
(`night`, `today`, `week`, `friends`) shared by all topics. The first
two topics of a student are drawn at random. After two different topics,
the next topic is the third one with probability 0.8 and otherwise drawn
at random; after a repeated topic, it is drawn at random. A first-order
model of the topic sequence cannot represent this rule.

## Usage

``` r
forum_posts
```

## Format

A data frame with 360 rows (one per post) and 5 columns:

- post:

  Character. The post identifier, `st01_t01` to `st30_t12`.

- student:

  Character. The author, `st01` to `st30`.

- turn:

  Integer. The position of the post among the author's posts, 1 to 12.

- topic:

  Character. The topic drawn by the simulation.

- text:

  Character. The text of the post, eight words.

## Source

Simulated with seed 11 by `data-raw/forum_posts.R`.

## See also

[`hg_sequences()`](https://mohsaqr.github.io/hypergraphs/reference/hg_sequences.md),
which turns a clustering of these posts into one sequence of topics per
student.

## Examples

``` r
head(forum_posts)
#>       post student turn     topic
#> 1 st01_t01    st01    1 astronomy
#> 2 st01_t02    st01    2 astronomy
#> 3 st01_t03    st01    3   cooking
#> 4 st01_t04    st01    4 gardening
#> 5 st01_t05    st01    5   cooking
#> 6 st01_t06    st01    6 astronomy
#>                                                           text
#> 1        stars nebula galaxy telescope moon nebula comet today
#> 2 stars stars galaxy galaxy stars astronomers astronomers week
#> 3           onions bake garlic oven garlic garlic bake friends
#> 4      watering tomatoes seeds prune soil prune tomatoes today
#> 5                   soup salt onions oven soup soup salt today
#> 6          orbit comet comet orbit comet telescope planet week
forum_hg <- text_hypergraph(forum_posts, column = "text", id = "post")
forum_hg
#> Text hypergraph: 360 documents, 34 words (documents as nodes, weight = n)
#> Hyperedges: 34 (words); sizes 48-104, median 63.5
#>       doc        word count weight
#>  st01_t01       comet     1      1
#>  st01_t01      galaxy     1      1
#>  st01_t01        moon     1      1
#>  st01_t01      nebula     2      2
#>  st01_t01       stars     1      1
#>  st01_t01   telescope     1      1
#>  st01_t01       today     1      1
#>  st01_t02 astronomers     2      2
#>  st01_t02      galaxy     2      2
#>  st01_t02       stars     3      3
#> ... 2229 more rows
```

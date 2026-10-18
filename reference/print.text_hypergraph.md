# Print a text hypergraph

Print a text hypergraph

## Usage

``` r
# S3 method for class 'text_hypergraph'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A
  [`text_hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/text_hypergraph.md)
  object.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
                        c = "stars and salt"))
print(hg, n = 3)
#> Text hypergraph: 3 documents, 4 words (documents as nodes, weight = n)
#> Hyperedges: 4 (words); sizes 2-3, median 2
#>  doc word count weight
#>    a  and     1      1
#>    a salt     1      1
#>    a soup     1      1
#> ... 6 more rows
```

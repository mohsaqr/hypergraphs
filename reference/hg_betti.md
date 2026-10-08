# Betti numbers of a simplicial complex

The ranks of the homology groups with rational coefficients: \\\beta_0\\
counts connected components, \\\beta_1\\ independent loops, \\\beta_2\\
voids, and so on.

## Usage

``` r
hg_betti(sc)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`simplicial()`](https://mohsaqr.github.io/hypergraphs/reference/simplicial.md).

## Value

A named integer vector `c(b0 = , b1 = , ...)`.

## Details

The coefficients are rational, not Z/2. The two agree on most complexes
but not on one with torsion: on the six-vertex real projective plane
(RP2) `hg_betti()` gives `(1, 0, 0)`, while
[`hg_homology()`](https://mohsaqr.github.io/hypergraphs/reference/hg_homology.md),
which reduces over Z/2 as the persistence literature does, ends its
Betti curve at `(1, 1, 1)`. Use
[`hg_homology()`](https://mohsaqr.github.io/hypergraphs/reference/hg_homology.md)
when Z/2 Betti numbers are meant.

## References

Hatcher, A. (2002). *Algebraic Topology*. Cambridge University Press.

## See also

[`hg_euler()`](https://mohsaqr.github.io/hypergraphs/reference/hg_euler.md),
[`hg_homology()`](https://mohsaqr.github.io/hypergraphs/reference/hg_homology.md)

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
hg_betti(simplicial(adj))
#> b0 b1 b2 
#>  1  0  0 
```

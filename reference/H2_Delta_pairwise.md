# Calculate pairwise heritability of differences between genotypes from model object

Instead of computing heritability on a "entry-mean" basis, this method
calculates heritability using "entry-differences". Entry here is
referring to the genotype, line or variety of interest. See reference
for origin and interpretation of `h2/H2_Delta_pairwise` and it's
variants

## Usage

``` r
H2_Delta_pairwise(model, target, type = c("BLUE", "BLUP"), options)
```

## Value

A `dspMatrix`

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## See also

[`H2_Delta_by_genotype()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_by_genotype.md),
[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md)

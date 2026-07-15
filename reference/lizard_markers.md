# Molecular markers data and genomic relatedness matrix of Australian lizards

Molecular markers data and genomic relatedness matrix of Australian
lizards

## Usage

``` r
lizard_markers

lizard_GRM
```

## Format

### `lizard_markers`

A data frame with 261 rows and 8438 columns:

- `gen` individual identifier, equivalent to `liz_id` in
  `lizard_phenotypes`

- 8437 genetic markers scored as 0, 2, 1 (see Details)

### `lizard_GRM`

A matrix array with 261 rows and 261 columns where each row/column
represents a genotype

## Source

<https://osf.io/hjkxd/?view_only=12a6b6010567474fac9fecd54472aa3d>

## Details

See `?snpReady::G.matrix` for how GRM was computed. The biallelic marker
\\M_iw\\ for the ith genotype and the \\w\\th marker with alleles
\\A_1\\ (i.e. the reference allele) and \\A_2\\ was coded as:

- 0 for \\A_1\\ \\A_1\\,

- 1 for \\A_1\\ \\A_2\\ and \\A_2\\ \\A_1\\

- 2 for \\A_2\\ \\A_2\\

## References

Kar, F., Nakagawa, S. and Noble, D.W.A. (2024) “Heritability and
developmental plasticity of growth in an oviparous lizard,” Heredity,
132(2), pp. 67–76. Available at:
https://doi.org/10.1038/s41437-023-00660-3.

# Calculate Standard heritability using variance parameters

Compute Standard heritability for genotype means using the variance
components of genotype and residuals.

## Usage

``` r
H2_Standard_parameters(vc_g, vc_e, n_r = 1)
```

## Arguments

- vc_g:

  Numeric. Genotype variance component

- vc_e:

  Numeric. Residuals variance component

- n_r:

  A numeric vector of size n_g, the number of genotype replicates.

## Value

Single numeric value

## Details

The equation for Standard heritability is as follows:

\$\$H^2\_{Standard} = \frac{\sigma^2_g}{\sigma^2_g +
\frac{1}{n_g}\sum\_{n_g}^{i=1} \sigma^2_p / n\_{gi}}\$\$ where:

- \\n_g\\ is the number of genotypes

- \\n\_{gi}\\ is the number of replicate for a given genotype i

- \\\sigma_g\\ is the variance attributed to genotype differences

- \\\sigma_p\\ is the variance attributed to phenotypic differences

## References

Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative
genetics (4th ed.). Longman.

## Examples

``` r
H2_Standard_parameters(vc_g = 0.25, vc_e = 0.8)
#> [1] 0.2380952
```

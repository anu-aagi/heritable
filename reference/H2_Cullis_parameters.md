# Calculate Cullis heritability using variance parameters

Compute the Cullis heritability for genotype means using the average
variance of pairwise differences of best linear unbiased predictors
(BLUPs).

## Usage

``` r
H2_Cullis_parameters(vd_BLUP_avg, vc_g)
```

## Arguments

- vd_BLUP_avg:

  Numeric. Average variance of pairwise differences among BLUPs

- vc_g:

  Numeric. Genotype variance component

## Value

Single numeric value

## Details

The equation for Cullis heritability is as follow

\$\$H^2\_{Cullis} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
ij}}{2\sigma^2_g}\$\$

where:

- \\PEV\\ is the prediction error variance matrix of the pairwise
  differences among BLUPS

- \\\sigma^2\\ is the variance attributed to differences between
  genotype

## References

Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of
early generation variety trials with correlated data. Journal of
Agricultural, Biological, and Environmental Statistics, 11(4), 381–393.
https://doi.org/10.1198/108571106X154443

## Examples

``` r
H2_Cullis_parameters(vd_BLUP_avg = 0.25, vc_g = 0.8)
#> [1] 0.84375
```

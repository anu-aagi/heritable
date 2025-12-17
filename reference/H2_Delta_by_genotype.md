# Calculate heritability of differences for a given genotype from model object

Instead of computing heritability on a "entry-mean" basis, this method
calculates heritability using "entry-differences". Entry here is
referring to the genotype, line or variety of interest. See reference
for origin and interpretation of `h2/H2_Delta_by_genotype` and it's
variants

## Usage

``` r
H2_Delta_by_genotype(model, target, type = c("BLUE", "BLUP"), options)
```

## Arguments

- model:

  Model object of class `lmerMod/merMod` or `asreml`

- target:

  The name of the random effect for which heritability is to be
  calculated.

- type:

  character, whether heritability is calculated using BLUEs or BLUPs

- options:

  NULL by default, for internal checking of model object before
  calculations

## Value

Numeric

Named list, with each element containing a named numeric vector

## Details

The heritability of differences for a given genotype is given by:

\$\$H^2\_{\Delta i.} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
i.}}{2\sigma^2_g}\$\$

where:

- \\PEV^{BLUP}\_{\overline\Delta i.}\\ is the arithmetic mean of the
  prediction error variance matrix for pairwise differences among BLUPs
  (or BLUEs if `method = "BLUE"`) for a given genotype

- \\\sigma^2\\ is the variance attributed to differences between
  genotype

See reference page 995 - 997 for full derivation of this heritability
measure and related variants

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## See also

[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md),
[`H2_Delta_pairwise()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_pairwise.md)

## Examples

``` r
# lme4 model
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Delta_by_genotype(lettuce_lme4, target = "gen", type = "BLUP")
#> $G1
#>        H2D_i
#> G1 0.8294971
#> 
#> $G10
#>         H2D_i
#> G10 0.8294971
#> 
#> $G11
#>         H2D_i
#> G11 0.8294971
#> 
#> $G12
#>         H2D_i
#> G12 0.8294971
#> 
#> $G13
#>         H2D_i
#> G13 0.8294971
#> 
#> $G14
#>         H2D_i
#> G14 0.8294971
#> 
#> $G15
#>         H2D_i
#> G15 0.8294971
#> 
#> $G16
#>         H2D_i
#> G16 0.8294971
#> 
#> $G17
#>         H2D_i
#> G17 0.8294971
#> 
#> $G18
#>         H2D_i
#> G18 0.8294971
#> 
#> $G19
#>         H2D_i
#> G19 0.8294971
#> 
#> $G2
#>        H2D_i
#> G2 0.8294971
#> 
#> $G20
#>         H2D_i
#> G20 0.8294971
#> 
#> $G21
#>         H2D_i
#> G21 0.8294971
#> 
#> $G22
#>         H2D_i
#> G22 0.8294971
#> 
#> $G23
#>         H2D_i
#> G23 0.8294971
#> 
#> $G24
#>         H2D_i
#> G24 0.8294971
#> 
#> $G25
#>         H2D_i
#> G25 0.8294971
#> 
#> $G26
#>         H2D_i
#> G26 0.8294971
#> 
#> $G27
#>         H2D_i
#> G27 0.8294971
#> 
#> $G28
#>         H2D_i
#> G28 0.8294971
#> 
#> $G29
#>         H2D_i
#> G29 0.8294971
#> 
#> $G3
#>        H2D_i
#> G3 0.8294971
#> 
#> $G30
#>         H2D_i
#> G30 0.8294971
#> 
#> $G31
#>         H2D_i
#> G31 0.8294971
#> 
#> $G32
#>         H2D_i
#> G32 0.8294971
#> 
#> $G33
#>         H2D_i
#> G33 0.8294971
#> 
#> $G34
#>         H2D_i
#> G34 0.8294971
#> 
#> $G35
#>         H2D_i
#> G35 0.8294971
#> 
#> $G36
#>         H2D_i
#> G36 0.8294971
#> 
#> $G37
#>         H2D_i
#> G37 0.8294971
#> 
#> $G38
#>         H2D_i
#> G38 0.8294971
#> 
#> $G39
#>         H2D_i
#> G39 0.8294971
#> 
#> $G4
#>        H2D_i
#> G4 0.8294971
#> 
#> $G40
#>         H2D_i
#> G40 0.8294971
#> 
#> $G41
#>         H2D_i
#> G41 0.8294971
#> 
#> $G42
#>         H2D_i
#> G42 0.8294971
#> 
#> $G43
#>         H2D_i
#> G43 0.8294971
#> 
#> $G44
#>         H2D_i
#> G44 0.8294971
#> 
#> $G45
#>         H2D_i
#> G45 0.8294971
#> 
#> $G46
#>         H2D_i
#> G46 0.8294971
#> 
#> $G47
#>         H2D_i
#> G47 0.8294971
#> 
#> $G48
#>         H2D_i
#> G48 0.8294971
#> 
#> $G49
#>         H2D_i
#> G49 0.8294971
#> 
#> $G5
#>        H2D_i
#> G5 0.8294971
#> 
#> $G50
#>         H2D_i
#> G50 0.8294971
#> 
#> $G51
#>         H2D_i
#> G51 0.8294971
#> 
#> $G52
#>         H2D_i
#> G52 0.8294971
#> 
#> $G53
#>         H2D_i
#> G53 0.8294971
#> 
#> $G54
#>         H2D_i
#> G54 0.8294971
#> 
#> $G55
#>         H2D_i
#> G55 0.8294971
#> 
#> $G56
#>         H2D_i
#> G56 0.8294971
#> 
#> $G57
#>         H2D_i
#> G57 0.8294971
#> 
#> $G58
#>         H2D_i
#> G58 0.8294971
#> 
#> $G59
#>         H2D_i
#> G59 0.8294971
#> 
#> $G6
#>        H2D_i
#> G6 0.8294971
#> 
#> $G60
#>         H2D_i
#> G60 0.8294971
#> 
#> $G61
#>         H2D_i
#> G61 0.8294971
#> 
#> $G62
#>         H2D_i
#> G62 0.8294971
#> 
#> $G63
#>         H2D_i
#> G63 0.8294971
#> 
#> $G64
#>         H2D_i
#> G64 0.8294971
#> 
#> $G65
#>         H2D_i
#> G65 0.8294971
#> 
#> $G66
#>         H2D_i
#> G66 0.8294971
#> 
#> $G67
#>         H2D_i
#> G67 0.8294971
#> 
#> $G68
#>         H2D_i
#> G68 0.8294971
#> 
#> $G69
#>         H2D_i
#> G69 0.8294971
#> 
#> $G7
#>        H2D_i
#> G7 0.8294971
#> 
#> $G70
#>         H2D_i
#> G70 0.8294971
#> 
#> $G71
#>         H2D_i
#> G71 0.8294971
#> 
#> $G72
#>         H2D_i
#> G72 0.8294971
#> 
#> $G73
#>         H2D_i
#> G73 0.8294971
#> 
#> $G74
#>         H2D_i
#> G74 0.8294971
#> 
#> $G75
#>         H2D_i
#> G75 0.8294971
#> 
#> $G76
#>         H2D_i
#> G76 0.8294971
#> 
#> $G77
#>         H2D_i
#> G77 0.8294971
#> 
#> $G78
#>         H2D_i
#> G78 0.8294971
#> 
#> $G79
#>         H2D_i
#> G79 0.8294971
#> 
#> $G8
#>        H2D_i
#> G8 0.8294971
#> 
#> $G80
#>         H2D_i
#> G80 0.8294971
#> 
#> $G81
#>         H2D_i
#> G81 0.8294971
#> 
#> $G82
#>         H2D_i
#> G82 0.8294971
#> 
#> $G83
#>         H2D_i
#> G83 0.8294971
#> 
#> $G84
#>         H2D_i
#> G84 0.8294971
#> 
#> $G85
#>         H2D_i
#> G85 0.8294971
#> 
#> $G86
#>         H2D_i
#> G86 0.8294971
#> 
#> $G87
#>         H2D_i
#> G87 0.8294971
#> 
#> $G88
#>         H2D_i
#> G88 0.8294971
#> 
#> $G89
#>         H2D_i
#> G89 0.8294971
#> 
#> $G9
#>        H2D_i
#> G9 0.8294971
#> 

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Delta_by_genotype(lettuce_asreml, target = "gen", type = "BLUP")
} # }
```

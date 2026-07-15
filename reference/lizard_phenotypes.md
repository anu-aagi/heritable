# Phenotypic data of Australian lizards

Phenotypic data for Australian lizards (*L. delicata*). Lizard eggs were
incubated at two different incubation treatment temperatures. Lizards
were weighed to investigate differences in growth between treatments.

## Usage

``` r
lizard_phenotypes
```

## Format

### `lizard_phenotypes`

A data frame with 3928 rows and 9 columns:

- `liz_id` individual identifier, equivalent to `gen` in
  `lizard_markers` 261 individuals

- `treatment` temperature (degrees celcius) at which lizard eggs were
  incubated

- `dam_id` dam identifier

- `sire_id` sire identifer

- `mass` mass of lizard in

- `ln_mass` natural log `mass`

- `days_since_hatch` days since hatching (age)

- `z_days_since_hatch` scaled and centered `days_since_hatch` see
  ?scale()

- `z_days_since_hatch_i2` squared `z_days_since_hatch` see Kar et al
  2024 to see details on their analysis

## Source

<https://osf.io/hjkxd/?view_only=12a6b6010567474fac9fecd54472aa3d>

## References

Kar, F., Nakagawa, S. and Noble, D.W.A. (2024) “Heritability and
developmental plasticity of growth in an oviparous lizard,” Heredity,
132(2), pp. 67–76. Available at:
https://doi.org/10.1038/s41437-023-00660-3.

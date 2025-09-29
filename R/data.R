#' Phenotypic and genetic marker data of 89 lettuce varieties
#'
#' @description
#' 89 lettuce varieties tested at three environments, each laid out as a randomized complete block design.
#' The measured trait was resistance to downy mildew scored on a scale ranging from 0 to 5.
#'
#' @format ## `lettuce_phenotypes`
#' A data frame with 703 rows and 4 columns:
#'
#' - `loc` environment identifier
#' - `gen` genotype identifier
#' - `rep` replicate identifier
#' - `y` resistance to downy mildew scored on a scale ranging from 0 to 5
#'
#' @format ## `lettuce_markers`
#' A data frame with 89 rows and 301 columns:
#'
#' - `gen` genotype identifier
#' - 300 genetic markers scored as -1, 0, 1 (see Details)
#' @references Hadasch, S., Simko, I., Hayes, R.J., Ogutu, J.O. and Piepho, H.-P. (2016), Comparing the Predictive Abilities of Phenotypic and Marker-Assisted Selection Methods in a Biparental Lettuce Population. The Plant Genome, 9: plantgenome2015.03.0014. [https://doi.org/10.3835/plantgenome2015.03.0014](https://doi.org/10.3835/plantgenome2015.03.0014)
#' @source [https://figshare.com/articles/dataset/Lettuce_trial_phenotypic_and_marker_data_/8299493](https://figshare.com/articles/dataset/Lettuce_trial_phenotypic_and_marker_data_/8299493)
#' @rdname lettuce
"lettuce_phenotypes"

#' @details The varieties were genotyped with a total of 300 markers (i.e. 95 single nucleotide polymorphisms and 205 amplified fragment length polymorphism markers, see Hayes et al. (2014) so that a marker matrix M_(89×300) was provided.
#' The biallelic marker M_iw for the ith genotype and the wth marker with alleles A_1 (i.e. the reference allele) and A_2 was coded as:
#'
#' - 1 for A_1 A_1,
#' -  -1 for A_2 A_2
#' - 0 for A_1 A_2 and A_2 A_1.
#' @rdname lettuce
#' @references Hayes, R. J., Galeano, C. H., Luo, Y., Antonise, R., & Simko, I. (2014). Inheritance of Decay of Fresh-cut Lettuce in a Recombinant Inbred Line Population from ‘Salinas 88’ × ‘La Brillante’. Journal of the American Society for Horticultural Science, 139(4), 388–398. [https://doi.org/10.21273/JASHS.139.4.388](https://doi.org/10.21273/JASHS.139.4.388)
"lettuce_markers"

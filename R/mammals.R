#' @title Example dataset for popbayes package
#'
#' @description ...
#'
#' @docType data
#'
#' @usage data(mammals)
#'
#' @format A data frame with 69 rows and 11 variables:
#'   - location: protected area name
#'   - species: species name
#'   - year: year
#'   - count: species counts
#'   - lower_bound: CI95% lower bound
#'   - upper_bound: CI95% upper bound
#'   - stat_method: S(ample), T(otal) or G(uesstimate)
#'   - field_method: A(erial) or G(round)
#'   - cv: coefficient of variation of species count
#'   - sd: standard deviation of species count
#'   - var: variance of species count
#'
#' @source \url{https://www.github.com/frbcesab/popbayes}
"mammals"

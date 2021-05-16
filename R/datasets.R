#' Counts conversion data
#' 
#' @description
#' This dataset contains information to convert individual counts estimated from
#' a field method to a preferred field method. The field method can be `A` 
#' (aerial counts) or `G` (ground counts). This dataset can be used by user to
#' add missing species. See [format_data()] for further information.
#'   
#' @format A data frame with 15 rows (species) and the following columns:
#' \describe{
#'   \item{order}{the order of the species}
#'   \item{family}{the family of the species}
#'   \item{species}{the species binomial name}
#'   \item{english}{the species English name}
#'   \item{french}{the species French name}
#'   \item{category}{the detectability category of the species. One of 
#'     \code{MLB} for Medium-sized Light and Brown species (20-150kg),
#'     \code{LLB} for Large Light and Brown species (>150kg),
#'     \code{LD} for Large Dark (>150kg), \code{Elephant}, and \code{Giraffe}}
#'   \item{pref_field_method}{the preferred field method of the species. One of 
#'     \code{A} for Aerial counts, and \code{G} for Ground counts}
#'   \item{conversion_fact}{the conversion factor (corresponding to the 
#'     detectability category) used to convert counts to the preferred field 
#'     method}
#' }
#' 
#' @examples 
#' data(conversion_data)
#' conversion_data

"conversion_data"

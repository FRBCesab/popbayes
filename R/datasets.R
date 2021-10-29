#' Species information dataset
#' 
#' @description
#' This dataset contains information about 15 African mammal species.
#' It can be used in the function `format_data()` to convert individual counts 
#' estimated from a field method to a preferred field method. The field method 
#' can be `A` (aerial counts) or `G` (ground counts). See [format_data()] for 
#' further information. It also contains the maximum population growth rate 
#' (i.e. the maximum change in log population size).
#' 
#' User can take this dataset as a template to add information for missing 
#' species. Note that only `species`, `pref_field_method`, `conversion_A2G`, 
#' and `rmax` are required.
#'   
#' @format A `data.frame` with 15 rows (African mammals species) and the 
#' following variables:
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
#'   \item{conversion_A2G}{the conversion multiplicative factor (corresponding  
#'     to the detectability category) used to convert aerial to ground counts}
#'   \item{rmax}{the maximum population growth rate}
#' }
#' 
#' @examples 
#' data("species_info")
#' species_info

"species_info"


#' African mammals survey in the Garamba National Park
#' 
#' @description
#' This dataset contains individual counts of 10 African mammal species
#' in the Garamba National Park (Democratic Republic of the Congo) from 1976 to
#' 2017.
#'   
#' @format A `data.frame` with 141 rows (counts) and the following 8 variables:
#' \describe{
#'   \item{location}{the location of the survey (Garamba)}
#'   \item{species}{the binomial name of the species}
#'   \item{date}{the date of the survey}
#'   \item{stat_method}{the method used to estimate individuals counts. One of
#'     \code{T} (total counts), \code{G} (guesstimate), and \code{S} 
#'     (sampling counts)}
#'   \item{field_method}{the field method used to collect data. One of 
#'     \code{A} (aerial counts), and \code{G} (ground counts)}
#'   \item{count}{number of individuals}
#'   \item{lower_ci}{lower boundary of the 95% confidence interval of the counts
#'     (only for sampling counts)}
#'   \item{upper_ci}{upper boundary of the 95% confidence interval of the counts
#'     (only for sampling counts)}
#'   \item{pref_field_method}{the preferred field method of the species. One of 
#'     \code{A} for Aerial counts, and \code{G} for Ground counts}
#'   \item{conversion_A2G}{the conversion multiplicative factor (corresponding  
#'     to the detectability category) used to convert aerial to ground counts}
#'   \item{rmax}{the maximum population growth rate}
#' }
#' 
#' @examples 
#' data("garamba")
#' head(garamba, 30)

"garamba"

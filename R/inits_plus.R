#' @title Generate Initial Values of Relative Rate of Increase
#'
#' @description
#' This function extracts the counts (c) and the corresponding dates (t) from 
#' the list of objects (data) to be analyzed by the bugs model. It calculates 
#' the observed annual relative rate of increase (r_cand_obs) for each > 0 
#' interval as the ratio of the change in log count divided by the interval 
#' length. It then draws a random initial value of the candidate relative rate 
#' of increase (rcand) to be passed to the bugs model as initial value by 
#' drawing in a normal distribution centered on the observed rate (r_cand_obs) 
#' with a standard deviation of 1. For each r_cand_obs, nc random values are 
#' drawn where nc (second argument) is the number of chains.
#'
#' @param data A named list of the data objects passed to the bugs model
#' @param nc A number of MCMC chains
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @importFrom stats rnorm
#' @export
#'
#' @return 
#' A list of nc lists (one per chain). Each list has exactly one element rcand, 
#' a vector of random initial values of candidate relative rates of increase, 
#' one per > 0 interval. For internal use only.
#'
#' @examples
#'
#' # No example

inits_plus <- function(data, nc){


  ### Arguments checks ----
  
  if (missing(data)) {
    stop("Argument `data` cannot be NULL.")
  }
  
  if (missing(nc)) {
    stop("Argument `nc` cannot be NULL.")
  }
  
  if (is.null(data)) {
    stop("Argument `data` cannot be NULL.")
  }

  if (is.null(nc)) {
    stop("Argument `nc` cannot be NULL.")
  }

  if (sum(names(data) == "c") == 0) {
    stop("Element `c` (counts) is missing from the list `data`.")
  }

  if (sum(names(data) == "t") == 0) {
    stop("Element `t` (time) is missing from the list `data`.")
  }

  if (sum(names(data) == "k") == 0) {
    stop("Element `k` (number of points) is missing from the list `data`.")
  }

  if (!is.numeric(data$c)) {
    stop("Element `c` (counts) of the list `data` must be numeric.")
  }

  if (length(data$c) < 4) {
    stop("Series must have at least 4 points.")
  }

  if (sum(data$c < 0) > 0) {
    stop("Element `c` (counts) of the list `data` must be positive or zero.")
  }

  if (!is.numeric(data$t)) {
    stop("Element `t` (time) of the list `data` must be numeric.")
  }

  if (sum(data$t <= 0) > 0) {
    stop("Element `t` (time) of the list `data` must be strictly positive.")
  }

  nc <- nc[1]

  if (!is.numeric(nc)) {
    stop("`nc` (number of MCMC chains) must be numeric of length 1.")
  }

  if (nc <= 0) {
    stop("`nc` (number of MCMC chains) must be > 0.")
  }

  
  ### Generate initial values of r ----
  
  r_cand_obs <- log(data$c[2:data$k]) - log(data$c[1:(data$k - 1)])
  l_interval <- data$t[2:data$k] - data$t[1:(data$k - 1)]

  position <- l_interval != 0
  r_cand_obs[position] <- r_cand_obs[position] / l_interval[position]

  list_start <- list()
  for (i in 1:nc) {
    list_start[[i]] <- list(
      rcand = stats::rnorm(length(r_cand_obs), r_cand_obs, 1)
    )
  }

  return(list_start)
}

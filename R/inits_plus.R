#' @title Generate initial values of relative rate of increase (r)
#'
#' @description
#' ...
#'
#' @param data [list] A named list of the data objects passed to the bugs model.
#' @param nc [integer] Number of MCMC chains.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' For internal use only.
#'
#' @examples
#'
#' No example

inits_plus <- function(data, nc){


  if (is.null(data)) {

    stop("`Data` cannot be NULL.")

  }

  if (is.null(nc)) {

    stop("`nc` cannot be NULL.")

  }


  if (sum(which(names(data) == "c")) == 0) {

    stop("Missing `c` element (counts) in `data`.")

  }

  if (sum(which(names(data) == "t")) == 0) {

    stop("Missing `t` element (time) in `data`.")

  }

  if (sum(which(names(data) == "k")) == 0) {

    stop("Missing `k` element (number of points) in `data`.")

  }

  if (!is.numeric(data$c)) {

    stop("Counts must be numeric.")

  }

  if (length(data$c) < 4) {

    stop("Series must have at least 4 points.")

  }

  if (sum(which(data$c < 0)) > 0) {

    stop("Counts must be positive or zero.")

  }

  if (!is.numeric(data$t)) {

    stop("Time must be numeric.")

  }

  if (sum(which(data$t <= 0)) > 0) {

    stop("Time must be strickly positive (years).")

  }


  nc <- nc[1]

  if (!is.numeric(nc)) {

    stop("`nc` (number of MCMC chains) must be numeric.")

  }

  if (nc <= 0) {

    stop("`nc` (number of MCMC chains) must be > 0.")

  }

  r_cand_obs <- log(data$c[2:data$k]) - log(data$c[1:(data$k - 1)])
  l_interval <- data$t[2:data$k] - data$t[1:(data$k - 1)]

  r_cand_obs[linterval != 0] <- r_cand_obs[l_interval != 0] / l_interval[l_interval != 0]

  list_start <- list()

  for (i in 1:nc) {

    list_start[[i]] <- list(rcand = rnorm(length(r_cand_obs), r_cand_obs, 1))
  }

  return(liststart)
}
